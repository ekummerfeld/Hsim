package edu.cmu.tetrad.algcomparison.independence;

import edu.cmu.tetrad.algcomparison.utils.Parameters;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.pitt.csb.mgm.IndTestMultinomialLogisticRegressionWald;

import java.util.Collections;
import java.util.List;

/**
 * Wrapper for Fisher Z test.
 *
 * @author jdramsey
 */
public class MultinomialLogisticRegressionWald implements IndependenceWrapper {

    @Override
    public edu.cmu.tetrad.search.IndependenceTest getTest(DataSet dataSet, Parameters parameters) {
        return new IndTestMultinomialLogisticRegressionWald(
                dataSet,
                parameters.getDouble("alpha"),
                false);
    }

    @Override
    public String getDescription() {
        return "MultinomialLogisticRetressionWald test";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        return Collections.singletonList("alpha");
    }

}
