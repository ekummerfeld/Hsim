package edu.cmu.tetrad.algcomparison.algorithm.continuous.dag;

import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.utils.Parameters;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Graph;

import java.util.ArrayList;
import java.util.List;

/**
 * LiNGAM.
 *
 * @author jdramsey
 */
public class Lingam implements Algorithm {
    public Graph search(DataSet dataSet, Parameters parameters) {
        edu.cmu.tetrad.search.Lingam lingam = new edu.cmu.tetrad.search.Lingam();
        return lingam.search(dataSet);
    }

    @Override
    public Graph getComparisonGraph(Graph graph) {
        return graph;
    }

    public String getDescription() {
        return "LiNGAM (Linear Non-Gaussian Acyclic Model";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        return new ArrayList<>();
    }
}
