package com.rtg.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMValidationError;
import java.util.List;
import junit.framework.TestCase;

public class Htsjdk6ValidationContractTest extends TestCase {
    public void testInvalidRecordStillSkipped() throws Exception {
        new SkipInvalidRecordsIteratorTest().checkWarning(
            SkipInvalidRecordsIteratorTest.SAM_HEAD10 + SkipInvalidRecordsIteratorTest.SAM_REC_BAD4,
            0,
            "WARNING::INVALID_FLAG_NOT_PRIM_ALIGNMENT",
            "WARNING::INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT",
            "At data line 1");
    }

    public void testUnmappedMappingQualityBoundary() {
        for (int quality : new int[]{0, 255, 1, 254}) {
            SAMRecord record = new SAMRecord(new SAMFileHeader());
            record.setReadName("mapping-quality-boundary");
            record.setReadUnmappedFlag(true);
            record.setMappingQuality(quality);
            List<SAMValidationError> errors = record.isValid();
            boolean invalidQuality = errors != null && errors.stream().anyMatch(
                error -> error.getType() == SAMValidationError.Type.INVALID_UNMAPPED_MAPPING_QUALITY);
            assertEquals("MAPQ=" + quality, quality != 0 && quality != 255, invalidQuality);
        }
    }
}
