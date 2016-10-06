package edu.cmu.tetrad.simulation;

import edu.cmu.tetrad.data.CovarianceMatrixOnTheFly;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Dag;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fgs;
import edu.cmu.tetrad.search.PatternToDag;
import edu.cmu.tetrad.search.SemBicScore;
import edu.cmu.tetrad.sem.*;

import java.io.PrintWriter;
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
public class SimEvalContinuous {
    public static void main(String... args){
        int verbosity = 2;
        //specify the number of iterations per data set type
        int iterations = 10;

        //specify the space of data params to span:
        List<Integer> numVars = Arrays.asList(15);
        List<Double> edgesPerNode = Arrays.asList(2.0);
        List<Integer> numCases = Arrays.asList(1000);

        //specify the space of predictor params to span:
        List<Integer> resimSize = Arrays.asList(1,3);
        List<Integer> hsimRepeat = Arrays.asList(40);
        List<Integer> fsimRepeat = Arrays.asList(40);

        String nl = System.lineSeparator();
        String output = "Simulation study output comparing Fsim and Hsim on predicting graph discovery accuracy"+nl;
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
                    //System.out.println(resimSize.size()+" "+hsimRepeat.size());
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
                        if (verbosity>1) System.out.println("iteration "+i);
                        Graph odag = GraphUtils.randomGraphRandomForwardEdges(varslist, 0, numEdges, 30, 15, 15, false, true);
                        //this focuses on nonlinear, nongaussian data
                        GeneralizedSemPm pm = HsimUtils.getNonlinearPm(odag,2);
                        GeneralizedSemIm im = new GeneralizedSemIm(pm);
                        //oData is the original data set, and odag is the original dag.
                        DataSet oData = im.simulateData(cases,false);
                        //now I have: a collection of data sets and graphs (maybe make an object for these pairs?)
                        //****optional: store the data and graphs somewhere separate
                        //**run FGS on data set, compare to graph. this makes the baseline for prediction.
                        // also used for fsim/hsim.
                        ICovarianceMatrix newcov = new CovarianceMatrixOnTheFly(oData);
                        SemBicScore oscore = new SemBicScore(newcov);
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
                                //then fit an IM to this dag and the data. GeneralizedSemEstimator seems to bug out
                                //GeneralizedSemPm simSemPm = new GeneralizedSemPm(fgsdag2);
                                //GeneralizedSemEstimator gsemEstimator = new GeneralizedSemEstimator();
                                //GeneralizedSemIm fittedIM = gsemEstimator.estimate(simSemPm, oData);

                                SemPm simSemPm = new SemPm(fgsdag2);
                                //BayesPm simBayesPm = new BayesPm(fgsdag2, bayesPm);
                                SemEstimator simSemEstimator = new SemEstimator(oData,simSemPm);
                                SemIm fittedIM = simSemEstimator.estimate();

                                DataSet simData = fittedIM.simulateData(cases, false);
                                //after making the full resim data (simData), run FGS on that
                                ICovarianceMatrix simcov = new CovarianceMatrixOnTheFly(simData);
                                SemBicScore simscore = new SemBicScore(simcov);
                                Fgs simfgs = new Fgs(simscore);
                                simfgs.setVerbose(false);
                                simfgs.setNumPatternsToStore(0);
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
                                HsimRepeatAC study = new HsimRepeatAC(oData);
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
                    String[][] latexTableArray = new String[resimSize.size()*hsimRepeat.size()+fsimRepeat.size()][5];
                    for (int j=0;j<fMSE.length;j++){
                        fMSE[j]=new PRAOerrors(fsimErrsByPars[j],"MSE for Fsim at vars="+vars+" edgeratio="+edgeratio+
                                " cases="+cases+" frepeat="+fsimRepeat.get(j));
                        if(verbosity>0){System.out.println(fMSE[j].allToString());}
                        output=output+fMSE[j].allToString()+nl;
                        latexTableArray[j]= prelimToPRAOtable(fMSE[j]);
                    }
                    for (int j=0;j<hMSE.length;j++){
                        for (int k=0;k<hMSE[j].length;k++){
                            hMSE[j][k]=new PRAOerrors(hsimErrsByPars[j][k],"MSE for Hsim at vars="+vars+" edgeratio="+edgeratio+
                                    " cases="+cases+" rsize="+resimSize.get(j)+" repeat="+hsimRepeat.get(k)+" iterations="+iterations);
                            if(verbosity>0){System.out.println(hMSE[j][k].allToString());}
                            output=output+hMSE[j][k].allToString()+nl;
                            latexTableArray[fsimRepeat.size()+j*hMSE[j].length+k]=prelimToPRAOtable(hMSE[j][k]);
                        }
                    }
                    //record all the params, the base error values, and the fsim/hsim mean squared errors
                    String latexTable = HsimUtils.makeLatexTable(latexTableArray);
                    try {
                        PrintWriter writer = new PrintWriter("latexTable.txt", "UTF-8");
                        writer.println(latexTable);
                        writer.close();
                    }
                    catch(Exception IOException){
                        IOException.printStackTrace();
                    }
                }
            }
        }
        try {
            PrintWriter writer = new PrintWriter("HvsF-SimulationEvaluation.txt", "UTF-8");
            writer.println(output);
            writer.close();
        }
        catch(Exception IOException){
            IOException.printStackTrace();
        }
    }
    //******************Private Methods************************//
    private static String[] prelimToPRAOtable(PRAOerrors input){
        String[] output = new String[5];
        double[] values = input.toArray();
        String[] vStrings = HsimUtils.formatErrorsArray(values,"%7.4e");
        output[0]=input.getName();
        for (int i=1;i<output.length;i++){
            output[i]=vStrings[i-1];
        }
        return output;
    }
}
