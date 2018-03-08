package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.python.StreamingPythonScriptExecutor;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.runtime.AsynchronousStreamWriterService;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.StreamSupport;

/**
 * Annotate a VCF with scores from a Convolutional Neural Network (CNN).
 *
 * This tool streams variants and their reference context to a python program
 * which evaluates a pre-trained neural network on each variant.
 * The neural network performs convolutions over the reference sequence surrounding the variant
 * and combines those features with a multilayer perceptron on the variant annotations.
 *
 * 2D models convolve over aligned reads as well as the reference sequence, and variant annotations.
 * 2D models require a SAM/BAM file as input and for the --tensor-type argument to be set
 * to a tensor type which requires reads, as in the example below.
 *
 * Pre-trained 1D and 2D models are included in the distribution.
 * It is possible to train your own models by using the tools:
 * {@link CNNVariantWriteTensors} and {@link CNNVariantTrain}.
 *
 *
 * <h3>1D Model Example</h3>
 *
 * <pre>
 * gatk CNNVariantScore\
 *   -V vcf_to_annotate.vcf.gz \
 *   -R reference.fasta \
 *   -O annotated.vcf \
 *   --architecture src/main/resources/org/broadinstitute/hellbender/tools/walkers/vqsr/cnn_1d_annotations.json
 * </pre>
 *
 * <h3>2D Model Example</h3>
 *
 * <pre>
 * gatk CNNVariantScore\
 *   -I aligned_reads.bam \
 *   -V vcf_to_annotate.vcf.gz \
 *   -R reference.fasta \
 *   -O annotated.vcf \
 *   -inference-batch-size 2 \
 *   -transfer-batch-size 2 \
 *   -tensor-type read-tensor \
 *   --architecture src/test/resources/large/VQSR/tiny_2d_wgs_tf_model.json
 * </pre>
 *
 */
@DocumentedFeature
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = CNNVariantScore.USAGE_SUMMARY,
        oneLineSummary = CNNVariantScore.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)

public class CNNVariantScore extends VariantWalker {
    private final static String NL = String.format("%n");
    static final String USAGE_ONE_LINE_SUMMARY = "Apply a Convolutional Neural Net to filter annotated variants";
    static final String USAGE_SUMMARY = "Annotate a VCF with scores from a Convolutional Neural Network (CNN)." +
            "The CNN determines a Log Odds Score for each variant." +
            "Pre-trained models (1D or 2D) are specified via the architecture argument." +
            "1D models will look at the reference sequence and variant annotations." +
            "2D models look at aligned reads, reference sequence, and variant annotations." +
            "2D models require a BAM file as input as well as the tensor-type argument to be set.";

    private static final int CONTIG_INDEX = 0;
    private static final int POS_INDEX = 1;
    private static final int REF_INDEX = 2;
    private static final int ALT_INDEX = 3;
    private static final int KEY_INDEX = 4;
    private static final int FIFO_STRING_INITIAL_CAPACITY = 1024;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file")
    private String outputFile;

    @Argument(fullName = "architecture", shortName = "architecture", doc = "Neural Net architecture configuration json file", optional = false)
    private String architecture;

    @Argument(fullName = "tensor-type", shortName = "tensor-type", doc = "Name of the tensors to generate, reference for 1D reference tensors and read_tensor for 2D tensors.", optional = true)
    private TensorType tensorType = TensorType.reference;

    @Argument(fullName = "window-size", shortName = "window-size", doc = "Neural Net input window size", minValue = 0, optional = true)
    private int windowSize = 128;

    @Advanced
    @Argument(fullName = "inference-batch-size", shortName = "inference-batch-size", doc = "Size of batches for python to do inference on.", minValue = 1, maxValue = 4096, optional = true)
    private int inferenceBatchSize = 256;

    @Advanced
    @Argument(fullName = "transfer-batch-size", shortName = "transfer-batch-size", doc = "Size of data to queue for python streaming.", minValue = 1, maxValue = 8192, optional = true)
    private int transferBatchSize = 512;

    @Hidden
    @Argument(fullName = "enable-journal", shortName = "enable-journal", doc = "Enable streaming process journal.", optional = true)
    private boolean enableJournal = false;

    @Hidden
    @Argument(fullName = "keep-temp-file", shortName = "keep-temp-file", doc = "Keep the temporary file that python writes scores to.", optional = true)
    private boolean keepTempFile = false;

    // Create the Python executor. This doesn't actually start the Python process, but verifies that
    // the requestedPython executable exists and can be located.
    final StreamingPythonScriptExecutor pythonExecutor = new StreamingPythonScriptExecutor(true);

    private FileOutputStream fifoWriter;
    private AsynchronousStreamWriterService<String> asyncWriter = null;
    private List<String> batchList = new ArrayList<>(inferenceBatchSize);

    private int curBatchSize = 0;
    private int windowEnd = windowSize / 2;
    private int windowStart = (windowSize / 2) - 1;
    private boolean waitforBatchCompletion = false;
    private File scoreFile;

    private String scoreKey;

    @Override
    protected String[] customCommandLineValidation() {
        if (inferenceBatchSize > transferBatchSize) {
            return new String[]{"Inference batch size must be less than or equal to transfer batch size."};
        }

        return null;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        scoreKey = getScoreKeyAndCheckModelAndReadsHarmony();

        // Start the Python process, and get a FIFO from the executor to use to send data to Python. The lifetime
        // of the FIFO is managed by the executor; the FIFO will be destroyed when the executor is destroyed.
        pythonExecutor.start(Collections.emptyList(), enableJournal);
        final File fifoFile = pythonExecutor.getFIFOForWrite();

        // Open the FIFO for writing. Opening a FIFO for read or write will block until there is reader/writer
        // on the other end, so before we open it, send an ASYNCHRONOUS command (one that doesn't wait for a
        // response) to the Python process to open the FIFO for reading. The Python process will then block until
        // we open the FIFO. We can then call getAccumulatedOutput.
        pythonExecutor.sendAsynchronousCommand(String.format("fifoFile = open('%s', 'r')" + NL, fifoFile.getAbsolutePath()));
        try {
            fifoWriter = new FileOutputStream(fifoFile);
        } catch (IOException e) {
            throw new GATKException("Failure opening FIFO for writing", e);
        }

        pythonExecutor.getAccumulatedOutput();
        asyncWriter = pythonExecutor.getAsynchronousStreamWriterService(fifoWriter, AsynchronousStreamWriterService.stringSerializer);
        batchList = new ArrayList<>(transferBatchSize);

        // Also, ask Python to open our output file, where it will write the contents of everything it reads
        // from the FIFO. <code sendSynchronousCommand/>
        try {
            scoreFile = File.createTempFile(outputFile, ".temp");
            if (!keepTempFile) {
                scoreFile.deleteOnExit();
            } else {
                logger.info("Saving temp file from python:" + scoreFile.getAbsolutePath());
            }
            pythonExecutor.sendSynchronousCommand(String.format("tempFile = open('%s', 'w+')" + NL, scoreFile.getAbsolutePath()));
            pythonExecutor.sendSynchronousCommand("import vqsr_cnn" + NL);
            pythonExecutor.sendSynchronousCommand(String.format("args, model = vqsr_cnn.args_and_model_from_semantics('%s')", architecture) + NL);
            logger.info("Using key:" + scoreKey + " for CNN architecture:" + architecture);
        } catch (IOException e) {
            throw new GATKException("Error when creating temp file and initializing python executor.", e);
        }

    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        referenceContext.setWindow(windowStart, windowEnd);
        if (tensorType.isReadsRequired()) {
            transferReadsToPythonViaFifo(variant, readsContext, referenceContext);
        } else {
            transferToPythonViaFifo(variant, referenceContext);
        }
        sendBatchIfReady();
    }

    private void transferToPythonViaFifo(final VariantContext variant, final ReferenceContext referenceContext) {
        try {
            final String outDat = String.format("%s\t%s\t%s\t%s\n",
                    getVariantDataString(variant),
                    new String(Arrays.copyOfRange(referenceContext.getBases(), 0, windowSize), "UTF-8"),
                    getVariantInfoString(variant),
                    variant.isSNP() ? "SNP" : variant.isIndel() ? "INDEL" : "OTHER");
            batchList.add(outDat);
            curBatchSize++;
        } catch (UnsupportedEncodingException e) {
            throw new GATKException("Trying to make string from reference, but unsupported encoding UTF-8.", e);
        }

    }

    private void sendBatchIfReady() {
        if (curBatchSize == transferBatchSize) {
            if (waitforBatchCompletion == true) {
                // wait for the last batch to complete before we start a new one
                asyncWriter.waitForPreviousBatchCompletion(1, TimeUnit.MINUTES);
                waitforBatchCompletion = false;
                pythonExecutor.getAccumulatedOutput();
            }
            executePythonCommand();
            waitforBatchCompletion = true;
            curBatchSize = 0;
            batchList = new ArrayList<>(transferBatchSize);
        }
    }

    private void transferReadsToPythonViaFifo(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext) {
        StringBuilder sb = new StringBuilder(FIFO_STRING_INITIAL_CAPACITY);
        try {
            sb.append(String.format("%s\t%s\t%s\t%s\t",
                    getVariantDataString(variant),
                    new String(Arrays.copyOfRange(referenceContext.getBases(), 0, windowSize), "UTF-8"),
                    getVariantInfoString(variant),
                    variant.isSNP() ? "SNP" : variant.isIndel() ? "INDEL" : "OTHER"));
        } catch (UnsupportedEncodingException e) {
            throw new GATKException("Trying to make string from reference, but unsupported encoding UTF-8.", e);
        }
        Iterator<GATKRead> readIt = readsContext.iterator();
        if (!readIt.hasNext()) {
            logger.warn("No reads at contig:" + variant.getContig() + "site:" + String.valueOf(variant.getStart()));
        }
        while (readIt.hasNext()) {
            sb.append(GATKReadToString(readIt.next()));
        }
        sb.append("\n");
        batchList.add(sb.toString());
        curBatchSize++;
    }

    private String GATKReadToString(GATKRead read) {
        StringBuilder sb = new StringBuilder(FIFO_STRING_INITIAL_CAPACITY);
        sb.append(read.getBasesString() + "\t");
        sb.append(baseQualityBytesToString(read.getBaseQualities()) + "\t");
        sb.append(read.getCigar().toString() + "\t");
        sb.append(read.isReverseStrand() + "\t");
        sb.append((read.isPaired() ? read.mateIsReverseStrand() : "false") + "\t");
        sb.append(read.isFirstOfPair() + "\t");
        sb.append(read.getMappingQuality() + "\t");
        sb.append(Integer.toString(read.getUnclippedStart()) + "\t");
        return sb.toString();

    }

    private String baseQualityBytesToString(byte[] qualities) {
        String qualityString = "";
        for (int i = 0; i < qualities.length; i++) {
            qualityString += Integer.toString(qualities[i]) + ",";
        }
        return qualityString.substring(0, qualityString.length() - 1);
    }

    private String getVariantDataString(final VariantContext variant) {
        return String.format("%s\t%d\t%s\t%s",
                variant.getContig(),
                variant.getStart(),
                variant.getReference().getBaseString(),
                variant.getAlternateAlleles().toString()
        );

    }

    private String getVariantInfoString(final VariantContext variant) {
        // Create a string that will easily be parsed as a python dictionary
        String varInfo = "";
        for (final String attributeKey : variant.getAttributes().keySet()) {
            varInfo += attributeKey + "=" + variant.getAttribute(attributeKey).toString().replace(" ", "").replace("[", "").replace("]", "") + ";";
        }
        return varInfo;
    }

    @Override
    public Object onTraversalSuccess() {
        if (waitforBatchCompletion) {
            asyncWriter.waitForPreviousBatchCompletion(1, TimeUnit.MINUTES);
            pythonExecutor.getAccumulatedOutput();
        }
        if (curBatchSize > 0) {
            executePythonCommand();
            asyncWriter.waitForPreviousBatchCompletion(1, TimeUnit.MINUTES);
            pythonExecutor.getAccumulatedOutput();
        }

        pythonExecutor.sendSynchronousCommand("tempFile.close()" + NL);
        pythonExecutor.sendSynchronousCommand("fifoFile.close()" + NL);
        pythonExecutor.terminate();

        writeOutputVCFWithScores();

        return true;
    }

    private void executePythonCommand() {
        final String pythonCommand = String.format(
                "vqsr_cnn.score_and_write_batch(args, model, tempFile, fifoFile, %d, %d)", curBatchSize, inferenceBatchSize) + NL;
        pythonExecutor.sendAsynchronousCommand(pythonCommand);
        asyncWriter.startAsynchronousBatchWrite(batchList);
    }


    private void writeOutputVCFWithScores() {
        try (final Scanner scoreScan = new Scanner(scoreFile);
             final VariantContextWriter vcfWriter = createVCFWriter(new File(outputFile))) {
            scoreScan.useDelimiter("\\n");
            writeVCFHeader(vcfWriter);
            final VariantFilter variantfilter = makeVariantFilter();

            // Annotate each variant in the input stream, as in variantWalkerBase.traverse()
            StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                    .filter(variantfilter)
                    .forEach(variant -> {
                        String sv = scoreScan.nextLine();
                        String[] scoredVariant = sv.split("\\t");
                        if (variant.getContig().equals(scoredVariant[CONTIG_INDEX])
                                && Integer.toString(variant.getStart()).equals(scoredVariant[POS_INDEX])
                                && variant.getReference().getBaseString().equals(scoredVariant[REF_INDEX])
                                && variant.getAlternateAlleles().toString().equals(scoredVariant[ALT_INDEX])) {
                            final VariantContextBuilder builder = new VariantContextBuilder(variant);
                            builder.attribute(scoreKey, scoredVariant[KEY_INDEX]);
                            vcfWriter.add(builder.make());
                        } else {
                            String errorMsg = "Score file out of sync with original VCF. Score file has:" + sv;
                            errorMsg += "\n But VCF has:" + variant.toStringWithoutGenotypes();
                            throw new GATKException(errorMsg);
                        }
                    });

        } catch (IOException e) {
            throw new GATKException("Error when trying to write annotated VCF.", e);
        }

    }

    private void writeVCFHeader(VariantContextWriter vcfWriter) {
        // setup the header fields
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();
        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
        hInfo.add(GATKVCFHeaderLines.getInfoLine(scoreKey));
        VariantRecalibrationUtils.addVQSRStandardHeaderLines(hInfo);
        final TreeSet<String> samples = new TreeSet<>();
        samples.addAll(inputHeader.getGenotypeSamples());
        hInfo.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    private String getScoreKeyAndCheckModelAndReadsHarmony() {
        if (tensorType.isReadsRequired() && this.hasReads()) {
            return GATKVCFConstants.CNN_2D_KEY;
        } else if (!tensorType.isReadsRequired()) {
            return GATKVCFConstants.CNN_1D_KEY;
        } else {
            throw new GATKException("2D Models require a SAM/BAM file specified via -I (-input) argument.");
        }
    }
}
