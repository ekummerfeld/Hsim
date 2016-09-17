package edu.cmu.tetrad.simulation;

import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Dag;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemImInitializationParams;
import edu.cmu.tetrad.sem.SemPm;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Clas for testing if HsimContinuous is working. Right now this basically just checks if
 * it can run without error on a basic data set.
 *
 * Created by Erich on 9/7/2016.
 */
public class TestHsimContinuous {
    public static void main(String... Args){
        int numEdges = 20;
        int cases = 50;
        int vars = 20;
        List<Node> varslist = new ArrayList<>();
        for (int i = 0; i < vars; i++) {
            varslist.add(new ContinuousVariable("X" + i));
        }
        Graph odag = GraphUtils.randomGraphRandomForwardEdges(varslist, 0, numEdges, 30, 15, 15, false, true);
        System.out.println(odag);
        SemIm oIM = new SemIm(new SemPm(odag),new SemImInitializationParams());
        DataSet oData = oIM.simulateDataRecursive(cases,false);
        System.out.println(oData);

        String[] resimNodeNames = {"X17","X18"};
        Set<Node> simnodes = new HashSet<Node>();
        for( int i = 0; i < resimNodeNames.length; i++) {
            Node thisNode = odag.getNode(resimNodeNames[i]);
            simnodes.add(thisNode);
        }
        Dag inputdag = new Dag(odag);
        System.out.println(inputdag);
        HsimContinuous hsim = new HsimContinuous(inputdag,simnodes,oData);
        DataSet hsimdata = hsim.hybridsimulate();
        System.out.println(hsimdata);
    }
}
