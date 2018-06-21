package org.broadinstitute.hellbender.engine;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.*;
import java.nio.file.spi.FileSystemProvider;
import java.util.HashMap;

/**
 * Default implementation for HtsURI.
 *
 * TODO: for simplicity, this currently contains some (cloud) functionality and tests that should ultimately
 * TODO: be moved into GatkURI, but its easier for now to see everything together
 * TODO: should we use the normalized uri string for hash/compare purposes, or the user-supplied string ?
 * TODO: special handling for authority/credentials, query params, fragment?
 *
 * TODO: fix FeatureInput.makeIntoAbsolutePath - makes a File Path out of a gendb path. amongst other things
 */
public class HtsSimpleURI implements HtsURI {

    private final String uriString;
    private URI uri;
    private Path cachedPath;
    private String reasonString;

    // There are 3 levels of validation:
    //
    // 1) HtsURI constructor - syntactically valid URI
    // 2) isNio - syntactically valid URI for which there is an installed NIO provider that matches the URI scheme
    // 3) isPath - syntactically valid URI that can be resolved to a Path by the provider
    //
    // <scheme>:<scheme-specific-part>
    // <scheme>://<authority><path>?<query>
    // absoluteURI   = scheme ":" ( hier_part | opaque_part )
    //      hier_part     = ( net_path | abs_path ) [ "?" query ]
    //      net_path      = "//" authority [ abs_path ]
    //      abs_path      = "/"  path_segments
    //
    // A URI is absolute if, and only if, it has a scheme component.
    // A URI is opaque if, and only if, it is absolute (has a scheme) and its
    // scheme-specific part does not begin with a slash character ('/')
    //
    //  A relative reference that does not begin with a scheme name or a
    //  slash character is termed a relative-path reference.
    //
    //  A relative reference beginning with a single slash character is
    //  termed an absolute-path reference, as defined by <abs_path> in
    //  Section 3.
    //
    // URI that do not make use of the slash "/" character for separating
    // hierarchical components are considered opaque by the generic URI
    // parser.
    //
    public HtsSimpleURI(final String uriString) {
        Utils.nonNull(uriString);
        this.uriString = uriString;

        try {
            uri = new URI(uriString);
            if (!uri.isAbsolute()) { // || uri.getScheme().equals("file")) {
                // No scheme component is present. Assume this a file URI; try to get a Path and then get
                // the URI from the resulting Path to get a properly escaped URI.
                //
                // If the input URI already has a file: scheme, then we assume that is is already properly
                // escaped.
                //
                // NOTE: This case (no scheme) is the only case where we resolve the URI to a Path at construction time.
                try {
                    setCachedPath(Paths.get(uriString));
                    uri = getCachedPath().toUri();
                } catch (InvalidPathException e) {
                    throw new IllegalArgumentException(e.getMessage(), e);
                }
            }
        } catch (URISyntaxException e) {
            final String errorMessage = String.format("%s must be a valid URI. '%s'/'%s'", uriString, e.getMessage(), e.getReason());
            throw new IllegalArgumentException(errorMessage);
        }

        //displayComponents(uri);
    }

    /**
     * Converts the URI to a {@link Path} object. If the filesystem cannot be found in the usual way, then attempt
     * to load the filesystem provider using the thread context classloader. This is needed when the filesystem
     * provider is loaded using a URL classloader (e.g. in spark-submit).
     *
     * Also makes an attempt to interpret the argument as a file name if it's not a URI.
     *
     * @return the resulting {@code Path}
     * @throws UserException if an I/O error occurs when creating the file system
     */
    @Override
    public Path toPath() {
        if (getCachedPath() != null) {
            return getCachedPath();
        }
        Path tmpPath;
        if (CloudStorageFileSystem.URI_SCHEME.equals(getURI().getScheme())) {
            tmpPath = BucketUtils.getPathOnGcs(getURI().toString());
        } else {
            tmpPath = Paths.get(getURI());
        }
        setCachedPath(tmpPath);
        return tmpPath;
    }

    @Override
    public URI getURI() {
        return uri;
    }

    @Override
    public String getURIString() {
        return uriString;
    }

    @Override
    public boolean isPath() {
        try {
            return getCachedPath() != null || toPath() != null;
        } catch (ProviderNotFoundException |
                FileSystemNotFoundException |
                IllegalArgumentException |
                UserException |
                AssertionError e) {
            // jimfs throws an AssertionError that wraps a URISyntaxException when trying to create path where
            // the scheme-specific part is missing or incorrect
            reasonString = e.getMessage();
            return false;
        }
    }

    @Override
    public boolean isNIO() {
        //if (!uri.isAbsolute()) {
        //    // !absolute == no scheme component. Assume this is a file on the default file system.
        //    return true;
        //}

        // try to find a provider
        for (FileSystemProvider provider: FileSystemProvider.installedProviders()) {
            if (provider.getScheme().equalsIgnoreCase(uri.getScheme())) {
                //setCachedPath(provider.getPath(uri));
                return true;
            }
        }
        return false;
    }

    @Override
    public String getToPathFailureReason() {
        if (reasonString == null) {
            try {
                toPath();
                return String.format("'%s' appears to be a valid Path", uriString);
            } catch (ProviderNotFoundException e) {
                return String.format("ProviderNotFoundException: %s", e.getMessage());
            } catch (FileSystemNotFoundException e) {
                return String.format("FileSystemNotFoundException: %s", e.getMessage());
            } catch (IllegalArgumentException e) {
                return String.format("IllegalArgumentException: %s", e.getMessage());
            } catch (UserException e) {
                return String.format("UserException: %s", e.getMessage());
            }
        }
        return reasonString;
    }

    // get the path associated with this URI
    protected Path getCachedPath() { return cachedPath; }

    protected void setCachedPath(Path path) {
        this.cachedPath = path;
    }

    @Override
    public String toString() {
        return uriString;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof HtsSimpleURI)) return false;

        HtsSimpleURI that = (HtsSimpleURI) o;

        if (!getURIString().equals(that.getURIString())) return false;
        if (!getURI().equals(that.getURI())) return false;
        return true;
    }

    @Override
    public int hashCode() {
        int result = getURIString().hashCode();
        result = 31 * result + getURI().hashCode();
        return result;
    }

    private void displayComponents(final URI uri) {
        System.out.println(
                String.format(
                        "Original:         %s\nURI:              %s\nASCII:            %s\nNormalized:       %s\nNormalized ASCII: %s\n",
                        uriString,
                        uri.toString(),
                        uri.toASCIIString(),
                        uri.normalize().toString(),
                        uri.normalize().toASCIIString()
                )
        );

        System.out.println(
                String.format(
                        "Scheme:         %s\nAuthority:      %s\nHost:           %s\n\n",
                        uri.getScheme(),
                        uri.getAuthority(),
                        uri.getHost()
                )
        );
    }

}
