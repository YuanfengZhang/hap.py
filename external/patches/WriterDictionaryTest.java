package com.rtg.sam;

import htsjdk.samtools.*;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import junit.framework.TestCase;

public class WriterDictionaryTest extends TestCase {
    public void testReferenceAndMateFollowOutputDictionary() throws Exception {
        SAMFileHeader input = new SAMFileHeader();
        input.addSequence(new SAMSequenceRecord("chr1", 100));
        input.addSequence(new SAMSequenceRecord("chr2", 100));
        SAMFileHeader output = new SAMFileHeader();
        output.addSequence(new SAMSequenceRecord("chr2", 100));
        output.addSequence(new SAMSequenceRecord("chr1", 100));
        SAMRecord record = new SAMRecord(input);
        record.setReadName("dictionary-remap");
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);
        record.setReferenceName("chr1");
        record.setAlignmentStart(2);
        record.setMateReferenceName("chr2");
        record.setMateAlignmentStart(3);
        record.setCigarString("1M");
        record.setReadString("A");
        record.setBaseQualityString("I");
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(output, true, bytes)) {
            writer.addAlignment(record);
        }
        try (SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(SamInputResource.of(new ByteArrayInputStream(bytes.toByteArray())))) {
            SAMRecord decoded = reader.iterator().next();
            assertEquals("chr1", decoded.getReferenceName());
            assertEquals("chr2", decoded.getMateReferenceName());
            assertEquals(2, decoded.getAlignmentStart());
            assertEquals(3, decoded.getMateAlignmentStart());
        }
    }
}
