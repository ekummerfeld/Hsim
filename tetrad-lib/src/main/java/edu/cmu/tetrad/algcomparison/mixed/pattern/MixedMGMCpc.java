package edu.cmu.tetrad.algcomparison.mixed.pattern;

import edu.cmu.tetrad.algcomparison.Algorithm;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.*;
import edu.pitt.csb.mgm.MGM;

import java.util.Map;

/**
 * Created by jdramsey on 6/4/16.
 */
public class MixedMGMCpc implements Algorithm {
    public Graph search(DataSet ds, Map<String, Number> parameters) {
        MGM m = new MGM(ds, new double[]{
                parameters.get("mgmParam1").doubleValue(),
                parameters.get("mgmParam2").doubleValue(),
                parameters.get("mgmParam3").doubleValue()
        });
        Graph gm = m.search();
        IndependenceTest indTest = new IndTestMixedLrt(ds, parameters.get("alpha").doubleValue());
        Cpc pcs = new Cpc(indTest);
        pcs.setDepth(-1);
        pcs.setInitialGraph(gm);
        pcs.setVerbose(false);
        return pcs.search();
    }

    public Graph getComparisonGraph(Graph dag) {
        return SearchGraphUtils.patternForDag(dag);
    }

    public String getDescription() {
        return "MGM-CPC: uses the output of MGM as an intial graph " +
                "for CPC-Stable, using the Mixed LRT test.";
    }
}