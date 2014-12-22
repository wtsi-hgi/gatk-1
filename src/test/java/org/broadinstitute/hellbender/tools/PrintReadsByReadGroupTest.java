/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.lang.reflect.Modifier;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.Callable;
import java.util.stream.Stream;

public class PrintReadsByReadGroupTest {
    private static final String TEST_DATA_PREFIX =
            "src/test/resources/org/broadinstitute/hellbender/tools/print_reads_by_read_group.";
    private static final File REFERENCE_SEQUENCE = new File(TEST_DATA_PREFIX + "fasta");

    @DataProvider(name = "printReadsByReadGroupData", parallel = true)
    public Object[][] getPrintReadsByReadGroupData() {
        return getSamReaderTypes()
                .filter(t -> t != SamReader.Type.CRAM_TYPE) // https://github.com/samtools/htsjdk/issues/148
                .map(t -> new Object[]{t, t == SamReader.Type.CRAM_TYPE})
                .toArray(Object[][]::new);
    }

    @Test(dataProvider = "printReadsByReadGroupData")
    public void testPrintReadsByReadGroup(final SamReader.Type type, final boolean useReference) throws Exception {
        final String fileExtension = type.fileExtension();
        final PrintReadsByReadGroup printReadsByReadGroup = new PrintReadsByReadGroup();
        printReadsByReadGroup.INPUT = new File(TEST_DATA_PREFIX + fileExtension);
        printReadsByReadGroup.OUTPUT_DIRECTORY = Files.createTempDirectory("printReadsByReadGroupData.");
        if (useReference)
            printReadsByReadGroup.REFERENCE_SEQUENCE = REFERENCE_SEQUENCE;
        printReadsByReadGroup.printReadsByReadGroup();

        Assert.assertEquals(
                getReadCounts(printReadsByReadGroup.OUTPUT_DIRECTORY, "Momma.0", fileExtension),
                17,
                "expected read group 0 count for " + fileExtension);
        Assert.assertEquals(
                getReadCounts(printReadsByReadGroup.OUTPUT_DIRECTORY, "Poppa.1", fileExtension),
                2,
                "expected read group 1 count for " + fileExtension);
    }

    private int getReadCounts(final Path tempDirectory, final String readGroupInfo, final String fileExtension) {
        final File path = tempDirectory.resolve(readGroupInfo + "." + fileExtension).toFile();
        int count = 0;
        IOUtil.assertFileIsReadable(path);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(path);
        for (@SuppressWarnings("unused") final SAMRecord rec : in) {
            count++;
        }
        CloserUtil.close(in);
        return count;
    }

    private static Stream<SamReader.Type> getSamReaderTypes() {
        return Stream
                .of(SamReader.Type.class.getFields())
                .filter(f -> Modifier.isStatic(f.getModifiers()))
                .filter(f -> f.getType().isAssignableFrom(SamReader.Type.class))
                .map(f -> orNull(() -> f.get(null)))
                .filter(v -> v instanceof SamReader.Type)
                .map(v -> (SamReader.Type) v);
    }

    private static <V> V orNull(final Callable<V> callable) {
        try {
            return callable.call();
        } catch (Exception e) {
            return null;
        }
    }
}
