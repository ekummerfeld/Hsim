package edu.cmu.tetrad.simulation;

import edu.cmu.tetrad.bayes.*;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Dag;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.BDeuScore;
import edu.cmu.tetrad.search.Fgs;
import edu.cmu.tetrad.search.PatternToDag;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by Erich on 8/22/2016.
 * //OUTLINE:
 //generate data sets and graphs (one loop for each param)
 //long loop generating howevermany datasets
 //****optional: store the data and graphs somewhere
 //now I have: a collection of data sets and graphs (maybe make an object for these pairs?)
 //iterate through that collection of pairs:
 //**run FGS on data set, compare to graph. this makes the baseline for prediction. also used for fsim/hsim.
 //**then step 1: full resim. iterate through the combinations of estimator parameters (just repeat num)
 //****calculate the squared errors of prediction, update a running sum of those values
 //**then step 2: hybrid sim. iterate through combos of params (repeat num, resimsize)
 //****calculate the squared errors of prediction, update running sum
 //divide the running sums by the size of the collection of pairs.
 //record all the params, the base error values, and the fsim/hsim mean squared errors
 */
public class SimulationEvaluation {
    public static void main(String... args){
        int verbosity = 1;
        //specify the number of iterations per data set type
        int iterations = 200;

        //specify the space of data params to span:
        List<Integer> numVars = Arrays.asList(7,12);
        List<Double> edgesPerNode = Arrays.asList(1.0);//0.5,1.0,1.5);
        List<Integer> numCases = Arrays.asList(200,800);

        //specify the space of predictor params to span:
        List<Integer> resimSize = Arrays.asList(2,3,5);
        List<Integer> hsimRepeat = Arrays.asList(10);//1,10,30,100);
        List<Integer> fsimRepeat = Arrays.asList(10);//1,10,30,100);

        //OUTLINE:
        //generate data sets and graphs (one loop for each param)
        for (Integer vars : numVars) {
            //initialize the set of variables
            List<Node> varslist = new ArrayList<>();
            for (int i = 0; i < vars; i++) {
                varslist.add(new DiscreteVariable("X" + i));
            }
            for (Double edgeratio : edgesPerNode){
                //calculate the number of edges
                final int numEdges = (int) (vars * edgeratio);
                for (Integer cases : numCases){
                    //*********8these store the errors from fsim and hsim, under different parameters
                    List<PRAOerrors>[] fsimErrsByPars = new ArrayList[fsimRepeat.size()];
                    int whichFrepeat = 0;
                    for (int frepeat : fsimRepeat){
                        fsimErrsByPars[whichFrepeat] = new ArrayList<PRAOerrors>();
                        whichFrepeat++;
                    }
                    List<PRAOerrors>[][] hsimErrsByPars = new ArrayList[resimSize.size()][hsimRepeat.size()];
                    System.out.println(resimSize.size()+" "+hsimRepeat.size());
                    int whichHrepeat;
                    int whichrsize = 0;
                    for (int rsize : resimSize) {
                        //hsimErrsByPars[whichrsize]=new ArrayList[hsimRepeat.size()];
                        whichHrepeat=0;
                        for (int hrepeat : hsimRepeat) {
                            //System.out.println(whichrsize+" "+whichHrepeat);
                            hsimErrsByPars[whichrsize][whichHrepeat]= new ArrayList<PRAOerrors>();
                            whichHrepeat++;
                        }
                        whichrsize++;
                    }
                    //*************done initializing hsimErrsByPars and fsimErrsByPars
                    //long loop generating howevermany datasets
                    for (int i=0;i<iterations;i++) {
                        Graph odag = GraphUtils.randomGraphRandomForwardEdges(varslist, 0, numEdges, 30, 15, 15, false, true);
                        BayesPm bayesPm = new BayesPm(odag, 2, 2);
                        BayesIm bayesIm = new MlBayesIm(bayesPm, MlBayesIm.RANDOM);
                        //oData is the original data set, and odag is the original dag.
                        DataSet oData = bayesIm.simulateData(cases, false);
                        //now I have: a collection of data sets and graphs (maybe make an object for these pairs?)
                        //****optional: store the data and graphs somewhere separate
                        //**run FGS on data set, compare to graph. this makes the baseline for prediction.
                        // also used for fsim/hsim.
                        BDeuScore oscore = new BDeuScore(oData);
                        Fgs ofgs = new Fgs(oscore);
                        ofgs.setVerbose(false);
                        ofgs.setNumPatternsToStore(0);
                        Graph oFGSGraph = ofgs.search();//***********This is the original FGS output on the data
                        PRAOerrors oErrors = new PRAOerrors(HsimUtils.errorEval(oFGSGraph, odag),"target errors");
                        if (verbosity>2) {
                            System.out.println(oFGSGraph);
                            System.out.println(oErrors.allToString());
                        }
                        //**then step 1: full resim. iterate through the combinations of estimator parameters (just repeat num)
                        for (whichFrepeat=0;whichFrepeat<fsimRepeat.size();whichFrepeat++){
                            ArrayList<PRAOerrors> errorsList = new ArrayList<PRAOerrors>();
                            for (int r =0;r<fsimRepeat.get(whichFrepeat);r++){
                                PatternToDag pickdag = new PatternToDag(oFGSGraph);
                                Graph fgsDag = pickdag.patternToDagMeek();
                                Dag fgsdag2 = new Dag(fgsDag);
                                BayesPm simBayesPm = new BayesPm(fgsdag2, bayesPm);
                                DirichletBayesIm simIM = DirichletBayesIm.symmetricDirichletIm(simBayesPm, 1.0);
                                DirichletEstimator simEstimator = new DirichletEstimator();
                                DirichletBayesIm fittedIM = simEstimator.estimate(simIM, oData);
                                DataSet simData = fittedIM.simulateData(cases, false);
                                BDeuScore simscore = new BDeuScore(simData);
                                Fgs simfgs = new Fgs(simscore);
                                simfgs.setVerbose(false);
                                simfgs.setNumPatternsToStore(0);
                                //simfgs.setPenaltyDiscount(penaltyDiscount);
                                Graph simGraphOut = simfgs.search();
                                PRAOerrors simErrors = new PRAOerrors(HsimUtils.errorEval(simGraphOut, fgsdag2), "Fsim errors "+r);
                                errorsList.add(simErrors);
                            }
                            PRAOerrors avErrors = new PRAOerrors(errorsList,"Average errors for Fsim at repeat="+fsimRepeat.get(whichFrepeat));
                            if (verbosity>3) System.out.println(avErrors.allToString());
                            //****calculate the squared errors of prediction, store all these errors in a list
                            double FsimAR2 = (avErrors.getAdjRecall()-oErrors.getAdjRecall()) *
                                    (avErrors.getAdjRecall()-oErrors.getAdjRecall());
                            double FsimAP2 = (avErrors.getAdjPrecision()-oErrors.getAdjPrecision()) *
                                    (avErrors.getAdjPrecision()-oErrors.getAdjPrecision());
                            double FsimOR2 = (avErrors.getOrientRecall()-oErrors.getOrientRecall()) *
                                    (avErrors.getOrientRecall()-oErrors.getOrientRecall());
                            double FsimOP2 = (avErrors.getOrientPrecision()-oErrors.getOrientPrecision()) *
                                    (avErrors.getOrientPrecision()-oErrors.getOrientPrecision());
                            PRAOerrors Fsim2 = new PRAOerrors(new double[]{FsimAR2,FsimAP2,FsimOR2,FsimOP2},
                                    "squared errors for Fsim at repeat="+fsimRepeat.get(whichFrepeat));
                            //add the fsim squared errors to the appropriate list
                            fsimErrsByPars[whichFrepeat].add(Fsim2);
                        }
                        //**then step 2: hybrid sim. iterate through combos of params (repeat num, resimsize)
                        for (whichrsize=0;whichrsize<resimSize.size();whichrsize++){
                            for (whichHrepeat=0;whichHrepeat<hsimRepeat.size();whichHrepeat++){
                                HsimRepeatAutoRun study = new HsimRepeatAutoRun(oData);
                                PRAOerrors HsimErrors= new PRAOerrors(study.run(resimSize.get(whichrsize), hsimRepeat.get(whichHrepeat)),"Hsim errors" +
                                        "at rsize="+resimSize.get(whichrsize)+" repeat="+hsimRepeat.get(whichHrepeat));
                                //****calculate the squared errors of prediction
                                double HsimAR2=(HsimErrors.getAdjRecall()-oErrors.getAdjRecall()) *
                                        (HsimErrors.getAdjRecall()-oErrors.getAdjRecall());
                                double HsimAP2=(HsimErrors.getAdjPrecision()-oErrors.getAdjPrecision()) *
                                        (HsimErrors.getAdjPrecision()-oErrors.getAdjPrecision());
                                double HsimOR2=(HsimErrors.getOrientRecall()-oErrors.getOrientRecall()) *
                                        (HsimErrors.getOrientRecall()-oErrors.getOrientRecall());
                                double HsimOP2=(HsimErrors.getOrientPrecision()-oErrors.getOrientPrecision()) *
                                        (HsimErrors.getOrientPrecision()-oErrors.getOrientPrecision());
                                PRAOerrors Hsim2 = new PRAOerrors(new double[]{HsimAR2,HsimAP2,HsimOR2,HsimOP2},
                                        "squared errors for Hsim, rsize="+resimSize.get(whichrsize)+" repeat="+hsimRepeat.get(whichHrepeat));
                                hsimErrsByPars[whichrsize][whichHrepeat].add(Hsim2);
                            }
                        }
                    }
                    //Average the squared errors for each set of fsim/hsim params across all iterations
                    PRAOerrors[] fMSE = new PRAOerrors[fsimRepeat.size()];
                    PRAOerrors[][] hMSE = new PRAOerrors[resimSize.size()][hsimRepeat.size()];
                    for (int j=0;j<fMSE.length;j++){
                        fMSE[j]=new PRAOerrors(fsimErrsByPars[j],"MSE for Fsim at (param details here)");
                        if(verbosity>0){System.out.println(fMSE[j].allToString());}
                    }
                    for (int j=0;j<hMSE.length;j++){
                        for (int k=0;k<hMSE[j].length;k++){
                            hMSE[j][k]=new PRAOerrors(hsimErrsByPars[j][k],"MSE for Hsim at (params)");
                            if(verbosity>0){System.out.println(hMSE[j][k].allToString());}
                        }
                    }
                    //record all the params, the base error values, and the fsim/hsim mean squared errors
                }
            }
        }

    }
}
