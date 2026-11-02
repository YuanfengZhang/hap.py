import com.fulcrumgenomics.jlibdeflate.LibdeflateCompressor;
import com.fulcrumgenomics.jlibdeflate.LibdeflateDecompressor;
import java.util.Arrays;
import java.util.Random;

public class NativeRoundtrip {
    public static void main(String[] args) {
        Random random = new Random(169);
        int cases = 0;
        for (int size : new int[]{0, 1, 1024, 65536, 1048576}) {
            byte[] input = new byte[size];
            random.nextBytes(input);
            try (var compressor = new LibdeflateCompressor();
                 var decompressor = new LibdeflateDecompressor()) {
                if (!Arrays.equals(input, decompressor.deflateDecompress(compressor.deflateCompress(input), size)))
                    throw new AssertionError("DEFLATE " + size);
                if (!Arrays.equals(input, decompressor.zlibDecompress(compressor.zlibCompress(input), size)))
                    throw new AssertionError("zlib " + size);
                if (!Arrays.equals(input, decompressor.gzipDecompress(compressor.gzipCompress(input), size)))
                    throw new AssertionError("gzip " + size);
                cases += 3;
            }
        }
        System.out.println("PASS native roundtrips: " + cases);
    }
}
