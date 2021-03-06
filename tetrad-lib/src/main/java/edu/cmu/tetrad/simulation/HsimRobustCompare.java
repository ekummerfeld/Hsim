package edu.cmu.tetrad.simulation;

/**
 * Created by Erich on 6/19/2016.
 */

import edu.cmu.tetrad.bayes.*;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataWriter;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.BDeuScore;
import edu.cmu.tetrad.search.Fgs;
import edu.cmu.tetrad.search.PatternToDag;
import edu.cmu.tetrad.util.RandomUtil;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * generate data from random graph, generated from parameters.
 * calculate errors from FGS output for the data and graph
 * create resimulated data and hybrid resimulated data with various parameters
 * calculate errors of FGS on the resimulated and hsim data
 * compare errors across all data sets. which simulated data errors are closest to original?
 */
public class HsimRobustCompare {

    //*************Public Methods*****************8//
    public static List<double[]> run(int numVars,double edgesPerNode, int numCases,double penaltyDiscount,
                              int resimSize,int repeat,boolean verbose) {
    //public static void main(String[] args) {
        //first generate the data
        RandomUtil.getInstance().setSeed(1450184147770L);
        char delimiter = ',';//'\t';
        final int numEdges = (int) (numVars * edgesPerNode);

        List<Node> vars = new ArrayList<>();
        double[] oErrors=new double[5];
        double[] hsimErrors=new double[5];
        double[] simErrors=new double[5];
        List<double[]> output = new ArrayList<>();

            for (int i = 0; i < numVars; i++) {
                vars.add(new DiscreteVariable("X" + i));
            }

            Graph odag = GraphUtils.randomGraphRandomForwardEdges(vars, 0, numEdges, 30, 15, 15, false, true);

            BayesPm bayesPm = new BayesPm(odag, 2, 2);
            BayesIm bayesIm = new MlBayesIm(bayesPm, MlBayesIm.RANDOM);
            //oData is the original data set, and odag is the original dag.
            DataSet oData = bayesIm.simulateData(numCases, false);
            //System.out.println(oData);
            //System.out.println(odag);

            //then run FGS
            BDeuScore oscore = new BDeuScore(oData);
            Fgs ofgs = new Fgs(oscore);
            ofgs.setVerbose(false);
            ofgs.setNumPatternsToStore(0);
            ofgs.setPenaltyDiscount(penaltyDiscount);
            Graph oGraphOut = ofgs.search();
            if (verbose) System.out.println(oGraphOut);

            //calculate FGS errors
            //oErrors = new double[5];
            oErrors = HsimUtils.errorEval(oGraphOut, odag);
            if (verbose) System.out.println(oErrors[0] + " " + oErrors[1] + " " + oErrors[2] +
                    " " + oErrors[3] + " " + oErrors[4]);

            //create various simulated data sets

            ////let's do the full simulated data set first: a dag in the FGS pattern fit to the data set.
            PatternToDag pickdag = new PatternToDag(oGraphOut);
            Graph fgsDag = pickdag.patternToDagMeek();

            Dag fgsdag2 = new Dag(fgsDag);
            BayesPm simBayesPm = new BayesPm(fgsdag2, bayesPm);
            DirichletBayesIm simIM = DirichletBayesIm.symmetricDirichletIm(simBayesPm, 1.0);
            DirichletEstimator simEstimator = new DirichletEstimator();
            DirichletBayesIm fittedIM = simEstimator.estimate(simIM, oData);
            DataSet simData = fittedIM.simulateData(numCases, false);

            ////next let's do a schedule of small hsims
            HsimRepeatAutoRun study = new HsimRepeatAutoRun(oData);
            hsimErrors = study.run(resimSize, repeat);

            //calculate errors for all simulated output graphs
            ////full simulation errors first
            BDeuScore simscore = new BDeuScore(simData);
            Fgs simfgs = new Fgs(simscore);
            simfgs.setVerbose(false);
            simfgs.setNumPatternsToStore(0);
            simfgs.setPenaltyDiscount(penaltyDiscount);
            Graph simGraphOut = simfgs.search();
            //simErrors = new double[5];
            simErrors = HsimUtils.errorEval(simGraphOut, fgsdag2);
            //System.out.println("Full resim errors are: " + simErrors[0] + " " + simErrors[1] + " " + simErrors[2] + " " + simErrors[3] + " " + simErrors[4]);

            //compare errors. perhaps report differences between original and simulated errors.
            //first, let's just see what the errors are.
            if (verbose)System.out.println("Original erors are: " + oErrors[0] + " " + oErrors[1] +
                    " " + oErrors[2] + " " + oErrors[3] + " " + oErrors[4]);
            if(verbose)System.out.println("Full resim errors are: " + simErrors[0] + " " + simErrors[1] +
                    " " + simErrors[2] + " " + simErrors[3] + " " + simErrors[4]);
            if(verbose)System.out.println("HSim errors are: " + hsimErrors[0] + " " + hsimErrors[1] +
                    " " + hsimErrors[2] + " " + hsimErrors[3] + " " + hsimErrors[4]);

            //then, let's try to squeeze these numbers down into something more tractable.
            //double[] ErrorDifferenceDifferences;
            //ErrorDifferenceDifferences = new double[5];
            //ErrorDifferenceDifferences[0] = Math.abs(oErrors[0]-simErrors[0])-Math.abs(oErrors[0]-hsimErrors[0]);
            //ErrorDifferenceDifferences[1] = Math.abs(oErrors[1]-simErrors[1])-Math.abs(oErrors[1]-hsimErrors[1]);
            //ErrorDifferenceDifferences[2] = Math.abs(oErrors[2]-simErrors[2])-Math.abs(oErrors[2]-hsimErrors[2]);
            //ErrorDifferenceDifferences[3] = Math.abs(oErrors[3]-simErrors[3])-Math.abs(oErrors[3]-hsimErrors[3]);
            //ErrorDifferenceDifferences[4] = Math.abs(oErrors[4]-simErrors[4])-Math.abs(oErrors[4]-hsimErrors[4]);
            //System.out.println("resim error errors - hsim error errors: " + ErrorDifferenceDifferences[0] + " " + ErrorDifferenceDifferences[1] + " " + ErrorDifferenceDifferences[2] + " " + ErrorDifferenceDifferences[3] + " " + ErrorDifferenceDifferences[4]);

        output.add(oErrors);
        output.add(simErrors);
        output.add(hsimErrors);

        return output;
    }
}