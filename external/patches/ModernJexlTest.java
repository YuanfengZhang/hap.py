package com.rtg.sam;
import junit.framework.TestCase;
import htsjdk.variant.variantcontext.VariantContextUtils;
import org.apache.commons.jexl3.MapContext;
public class ModernJexlTest extends TestCase {
  public void testEvaluationAndMissingValues() {
    MapContext context = new MapContext(); context.set("QUAL", 30);
    assertEquals(Boolean.TRUE, VariantContextUtils.engine.get().createExpression("QUAL >= 20").evaluate(context));
    assertEquals(Boolean.FALSE, VariantContextUtils.engine.get().createExpression("QUAL < 20").evaluate(context));
    try { VariantContextUtils.engine.get().createExpression("MISSING > 1").evaluate(context); fail("missing variable accepted"); } catch (org.apache.commons.jexl3.JexlException expected) { }
  }
}
