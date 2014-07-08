//Andrew Peterman
//Random walk with shortest path (Dijkstra's) between starting node and all pathogens
//Updated 2014 for k-core
//RRWBiasedWalk2 (TBRRW)

package TBRRW;
import java.io.BufferedReader;
import java.io.FileReader;
import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.Algebra;
import java.util.*;
public class Main { 
    //Parameters
 
    static double a = 0.6; //restart probabilty
    static double lamda = .55; //cutoff for cluster expansion
    static int bias = 1; 
    static int biasSteps = 3; 
    static int biasDegree = 0; 
    static int k = 20; //max cluster size 
    
    static ArrayList<String> ConnectionList =  new ArrayList<String>();
    static ArrayList<String> proteinList =  new ArrayList<String>();
    static double[][] transition;
    static SparseDoubleMatrix2D Pt;
    static ArrayList<Integer> pathogenList =  new ArrayList<Integer>();
    static ArrayList<Integer> degree =  new ArrayList<Integer>();
    static ArrayList<Double> nodeCoreDensity =  new ArrayList<Double>();
    static ArrayList<Integer> nodeCoreSize =  new ArrayList<Integer>();
    static ArrayList<Integer> toBias =  new ArrayList<Integer>();
    static HashSet toB = new HashSet();
    static ArrayList<ArrayList<String>> bestPath =  new ArrayList<ArrayList<String>>();
    static ArrayList<ArrayList<Integer>> kCores = new ArrayList<ArrayList<Integer>>();
    static long startTime = System.nanoTime();    
    //Main program 
    public static void main(String[] args) { 
/*******************************Build Transition Matrix************************/
        
        long endTime; 
                 
        //Read in Graph File
        try {
            BufferedReader readbuffer;
            //readbuffer = new BufferedReader(new FileReader("C:\\Users\\apeterman\\Downloads\\_apstuff\\890\\yeastW.txt"));
            //readbuffer = new BufferedReader(new FileReader("I:\\_apstuff\\_Research\\HIV_H_MSVM_1s.txt"));
            //readbuffer = new BufferedReader(new FileReader("E:\\_apstuff\\_Research\\VIFW1.txt"));
            //readbuffer = new BufferedReader(new FileReader("F:\\_apstuff\\_Research\\CUL5W.txt"));
            readbuffer = new BufferedReader(new FileReader("C:\\Users\\Joe Schmoe\\Downloads\\899\\TESTW2.txt"));
//            readbuffer = new BufferedReader(new FileReader(args[0])); 
//            a = Double.valueOf(args[1]);
//            lamda = Double.valueOf(args[2]);
//            bias = Integer.valueOf(args[3]);
//            biasSteps = Integer.valueOf(args[4]);
//            biasDegree = Integer.valueOf(args[5]);
            
            String strRead;
            //Process File build MxM matrix - find distinct proteins, connected proteins
            while ((strRead=readbuffer.readLine())!=null){
                //Keep the original list to make matrix
                ConnectionList.add(strRead);
                String splitarray[] = strRead.split("\t");
                    //Make distinct list of proteins
                    String firstProtein = splitarray[0].toUpperCase();
                    if(!(proteinList.contains(firstProtein)))
                    {
                        proteinList.add(firstProtein);
                    }
                    else
                    { 
                    }
                    String secondProtein = splitarray[1].toUpperCase();
                    if(!(proteinList.contains(secondProtein)))
                    {
                        proteinList.add(secondProtein);
                    }
                    else
                    {
                    }
                    //System.out.println(proteinList.indexOf(firstProtein) + " " + proteinList.indexOf(secondProtein));
            }
            readbuffer.close();
        } catch (Exception ex) {
            System.out.println("Error:" + ConnectionList.size());
        }
        System.out.println("Number of Proteins: " + proteinList.size() + "\n");
        
        //Build pathogen list
        for(int i=0;i<proteinList.size();i++)
            if (proteinList.get(i).equals("CA") || proteinList.get(i).equals("GAG")
                || proteinList.get(i).equals("GP120") || proteinList.get(i).equals("GP160")
                || proteinList.get(i).equals("GP41") || proteinList.get(i).equals("IN")
                || proteinList.get(i).equals("MA") || proteinList.get(i).equals("NC")
                || proteinList.get(i).equals("NEF") || proteinList.get(i).equals("POL")
                || proteinList.get(i).equals("P6") || proteinList.get(i).equals("P1")
                || proteinList.get(i).equals("PR") || proteinList.get(i).equals("REV")
                || proteinList.get(i).equals("RT") || proteinList.get(i).equals("TAT")
                || proteinList.get(i).equals("VIF") || proteinList.get(i).equals("VPR")
                || proteinList.get(i).equals("VPU"))
            {
                pathogenList.add(i);
            }
        
       //Make array for network edges
       transition = new double [proteinList.size()][proteinList.size()];
       //copy transition to keep weightings for later
       //double transNormalized[][] = new double [proteinList.size()][proteinList.size()];
       //fill the array

       String edges;
       String splitarray[];
       System.out.println("Creating transition matrix...");
       for(int i=0;i<ConnectionList.size();i++)
       {
           edges = ConnectionList.get(i);
           splitarray = edges.split("\t");
           String P1 = splitarray[0].toUpperCase();
           String P2 = splitarray[1].toUpperCase();
           double C = Double.valueOf(splitarray[2]);
           int c1 = proteinList.indexOf(P1);
           int c2 = proteinList.indexOf(P2);
           transition[c1][c2] = C;
           transition[c2][c1] = C;
       }
       
       StoreDegreeDistribution(transition);
       
       for(int i=0;i<proteinList.size();i++)
       {
           nodeCoreSize.add(0);
           nodeCoreDensity.add(0.0);  
       }       
       
       CalcKcore(degree);
       

       
       Pt = new SparseDoubleMatrix2D(transition);
       
       double transitionCopy[][] = new double [proteinList.size()][proteinList.size()];

        //deep copy transition array to bias
        for (int t = 0; t < transition.length; t++)
            transitionCopy[t] = Arrays.copyOf(transition[t], transition[t].length);

       //normalize probabilities
       for(int b=0;b<transitionCopy.length;b++)
       {         
            double sum=0;
            for(int j=0;j<transitionCopy[b].length;j++)
            {
                sum = sum + transitionCopy[j][b];
            }
            for(int j=0;j<transitionCopy[b].length;j++)
            {
                transitionCopy[j][b] = transitionCopy[j][b]/sum;
            }     
       }
       
       SparseDoubleMatrix2D P = new SparseDoubleMatrix2D(transitionCopy);

        endTime = System.nanoTime();
        System.out.println("Elaspsed Time: " + ((endTime-startTime)/1000000000.0) + " sec\n");
       //Get shortest pathogen paths to all nodes     
        System.out.println("Finding path lengths to all pathogens...");
        for(int i=0;i<pathogenList.size();i++)
        {
            getPathogenPaths(pathogenList.get(i));
        }
       //RRW using edge weights 
       endTime = System.nanoTime(); 
       System.out.println("Elaspsed Time: " + ((endTime-startTime)/1000000000.0) + " sec\n");    

       RepeatedRandomWalks(proteinList,P); 
    }
    //Calculate Random walk from a node with restart probability
    public static DenseDoubleMatrix1D RandomWalkWithRestart(SparseDoubleMatrix2D P, int length, int n, double a){
        
        Algebra algebra = new Algebra();
        //Restart vector, all zeros except for protein
        double[] restartVector = new double [length];
        restartVector[n] = a;
        DenseDoubleMatrix1D s = new DenseDoubleMatrix1D(restartVector);
        //Stationary vector
        restartVector[n] = 1;
        DenseDoubleMatrix1D x = new DenseDoubleMatrix1D(restartVector);
        //Comparison vector for norm
        DenseDoubleMatrix1D x1 = new DenseDoubleMatrix1D(restartVector);
        double[] nonRestart = new double[length];
        double val = (1.0-a);
        for(int i=0; i < length;i++) {
            nonRestart[i] = val;
        }
        DenseDoubleMatrix1D vC = new DenseDoubleMatrix1D(nonRestart);
        //Check L2 norm to 9 decimal places having not changed or stop at 75 iterations
        int iterations;
        for(iterations=0;iterations<75;iterations++)
        {
            //System.out.println(x);
            //Save x to compare in next iteration
            x1.assign(x);
            //Calc steady-state probability x = a*s + (1-a)*Pt*x
            x = (DenseDoubleMatrix1D) x.assign(vC,cern.jet.math.Functions.mult);
            x = (DenseDoubleMatrix1D) P.zMult(x, null).assign(s, cern.jet.math.Functions.plus);
            if(Math.abs(algebra.norm2(x1) - algebra.norm2(x)) < 0.0000000001)
                break;
        }
        //System.out.print(x); 
        return x;
    }
    
    //Average the steady-state vectors of C for the RRW method
    public static DenseDoubleMatrix1D ClusterWalks(ArrayList<Integer> C, ArrayList<DenseDoubleMatrix1D> m){
        double[] s = new double [m.size()];
        double[] s1 = new double [m.size()];
        //DenseDoubleMatrix1D sum = new Matrix(s, 1);
        DenseDoubleMatrix1D sum = new DenseDoubleMatrix1D(s);
        DenseDoubleMatrix1D sum1 = new DenseDoubleMatrix1D(s1);
        //initialize a
        double[] thing = new double[m.size()];
        double val =  1.0/C.size();
        for(int i=0; i < m.size();i++) {
            thing[i] = val;
        }
        DenseDoubleMatrix1D vC = new DenseDoubleMatrix1D(thing);
        for(int i=0;i<C.size();i++)
        {
            sum1 = (DenseDoubleMatrix1D) m.get(C.get(i));
            sum = (DenseDoubleMatrix1D) sum.assign(sum1,cern.jet.math.Functions.plus);
        }
        sum = (DenseDoubleMatrix1D) sum.assign(vC,cern.jet.math.Functions.mult);
        return sum;
    }
            
    //Bias paths in temporary transition matrix for random walk from n to all pathogens
    private static void BiasTransitionMatrix(int node, SparseDoubleMatrix2D bias) {

        //Update paths to HIV proteins
        //With list of path lengths we need to bias the edges of the shortest paths to the pathogens        
        //For a node, update all pathogen pathway
                                        
        ArrayList<String> nodePaths = new ArrayList<String>();
        nodePaths = bestPath.get(node);
        toBias.clear();
        if(nodePaths.size() > 0)
        {
            for(int i=0;i<nodePaths.size();i++)
            { 
                String edges = nodePaths.get(i);
                //System.out.println(edges);
                String[] h = edges.split(",");
                
                //don't bias long paths longer than this parameter
                if(h.length >= biasSteps+1)
                    break;
                
                //System.out.print(bias);
                //for each node that is along the path copy columns from original matrix, bias and normalize
                for(int b=0;b<h.length;b++)
                {                       
                    int currentNode = Integer.parseInt(h[b]);
                    //int nextNode = Integer.parseInt(h[b+1]);
                    //System.out.println("Updating " + currentNode + " and " + nextNode);

                    //copy column from transition
                    if(!toBias.contains(currentNode))
                        bias.viewColumn(currentNode).assign(Pt.viewColumn(currentNode));
                    
                    if(!toBias.contains(currentNode))
                        toBias.add(currentNode);
                    //System.out.print(bias);
                    
                }
                //
                for(int b=0;b<h.length-1;b++)
                {                                           
                    int currentNode = Integer.parseInt(h[b]);
                    int nextNode = Integer.parseInt(h[b+1]);
                    
                    if(biasDegree > 0)
                    {
                        int d1 = degree.get(currentNode);
                        int d2 = degree.get(nextNode);
                        if (d2 < d1)
                            break;
                    }
                    
                    //version 1
                    //double biasValue = nodeCoreSize.get(node) * nodeCoreDensity.get(node);
                    //version 2
                    double biasValue = 1.0/(h.length-1) * nodeCoreSize.get(node) * nodeCoreDensity.get(node);
                    
                    //bias the two node edges    
                    bias.setQuick(currentNode, nextNode, biasValue);
                    //bias[currentNode][nextNode] = 1;
                    bias.setQuick(nextNode, currentNode, biasValue);
                    //bias[nextNode][currentNode] = 1;
                }
                 //System.out.print(bias);
                 
                 //Updated rows can get re-updated in the second pass if
                 //there are more than 1 protein. The biasing may be lost 
                 //if the second path overlaps. Better to do it at the very end??
                
                //System.out.print(bias);      
            } 
            
            for(int b=0;b<toBias.size();b++)
            {   
                //double sum=0;
                int currentNode = toBias.get(b);
                bias.viewColumn(currentNode).assign(cern.jet.math.Functions.div(bias.viewColumn(currentNode).zSum()));

            }
        }
       //System.out.print(bias); 
    }
    
    //Bias paths in temporary transition matrix for random walk from n to all pathogens
    private static void RestoreTransitionMatrix(int node, double [][] bias) {

        //Update paths to HIV proteins
        //With list of path lengths we need to bias the edges of the shortest paths to the pathogens        
        //For a node, update all pathogen pathway
        
        ArrayList<String> nodePaths = new ArrayList<String>();
        nodePaths = bestPath.get(node);
        if(nodePaths.size() > 0)
        {
            for(int i=0;i<nodePaths.size();i++)
            {
                String edges = nodePaths.get(i);
                String[] h = edges.split(",");

                for(int b=0;b<h.length-1;b++)
                {   
                    bias[Integer.parseInt(h[b])][Integer.parseInt(h[b+1])] = 1;
                    bias[Integer.parseInt(h[b+1])][Integer.parseInt(h[b])] = 1;

                    if(!toBias.contains(b))
                        toBias.add(b);
                        //toB.add(b);
                }
            }
            
            for(int b=0;b<bias.length;b++)
            {
                if(toBias.contains(b))
                {
                    double sum=0;
                    for(int j=0;j<bias[b].length;j++)
                    {
                        sum = sum + bias[j][b];
                    }
                    for(int j=0;j<bias[b].length;j++)
                    {
                        bias[j][b] = bias[j][b]/sum;
                    } 
                }    
            }
                  
        }
       
    }    
    
    
    private static void getPathogenPaths(int pathogen) {
        
        //update the best paths to each node from the pathogen        
        int [] best = new int [transition.length];
        boolean[] visited = new boolean[transition.length];
        int max = 10000; // Infinity equivalent.
        String path;
                
        //create new path for this pathogen
        if(bestPath.size() > 0)
        {
            for (int i = 0; i < transition.length; i++)
            {     
                path = Integer.toString(pathogen)+ ",";
                ArrayList<String> pathList = new ArrayList<String>();

                pathList = bestPath.get(i);
                pathList.add(path);
                bestPath.set(i,pathList);        
            }
         }
        else
        {
            for (int i = 0; i < transition.length; i++)
            {     
                path = Integer.toString(pathogen)+ ",";
                ArrayList<String> pathList = new ArrayList<String>();
                
                pathList.add(path);
                bestPath.add(pathList);                
            }
        }
         
        //initalize
        for (int i = 0; i < transition.length; i++)
        {
            best[i] = max;
            visited[i] = false;
        }
        
        best[pathogen] = 0;
        //System.out.print("Starting node " + n + "\n");
        for(int i = 0; i < transition.length; i++)
        {
            double min = max;
            int currentNode = 0;
            
            //for each of the nodes that's not visited, get the shortest path to the other nodes
            for(int j = 0; j < transition.length; j++)
            {
                if(!visited[j] && best[j] < min) 
                {                    
                    currentNode = j; 
                    min = best[j];
                }
            } 
            //mark the best node as the one visited - they are all the same for unweighted
            visited[currentNode] = true;
            
            //update all distances from currentNode if they are better 
            for(int j = 0; j < transition.length; j++)
            {
                if ((transition[currentNode][j] != 0) && (best[currentNode] + 1) < best[j])
                { 
                    best[j] = best[currentNode] + 1;//transition[currentNode][j]; 
                    ArrayList<String> newPath =  new ArrayList<String>();
                    String newP;
                    newPath = bestPath.get(j);
                    //System.out.println(currentNode);
                    newP = bestPath.get(currentNode).get(pathogenList.indexOf(pathogen)) + Integer.toString(j) + ',';
                    newPath.set(pathogenList.indexOf(pathogen),newP);
                    bestPath.set(j,newPath);
                }              
            }
        }        
    }
    
    //Macropol Method
    public static void RepeatedRandomWalks(ArrayList<String> proteinList, SparseDoubleMatrix2D P){
       
        ArrayList<DenseDoubleMatrix1D> m = new ArrayList<DenseDoubleMatrix1D>();
        long endTime = System.nanoTime();
        //steps RRW 1-4     
        
//        double transitionBias[][] = new double [proteinList.size()][proteinList.size()];
//        if (bias ==1)  
//        {
//            //deep copy transition array to bias        
//            for (int t = 0; t < transNormalized.length; t++)
//                transitionBias[t] = Arrays.copyOf(transNormalized[t], transNormalized[t].length);        
//        }
        System.out.println("Creating All Random Walk Vectors...");
        //int count=0;
        for(int i=0;i<proteinList.size();i++)
        {
            //bias the transition array for this RWR
            if (bias ==1)
            {
                //dont copy....keep list, bias and un-bias the couple
                                               
                //endTime = System.nanoTime(); 
                //System.out.println("Before Bias " + ((endTime-startTime)/1000000000.0) + " sec\n");
                //SparseDoubleMatrix2D P = new SparseDoubleMatrix2D(transNormalized);
                BiasTransitionMatrix(i,P); 
                //endTime = System.nanoTime();
                //System.out.println("Before Normalize " + ((endTime-startTime)/1000000000.0) + " sec\n");                

                //normalize the updated columns

                 
                //endTime = System.nanoTime();
                //System.out.println("Before Walk " + ((endTime-startTime)/1000000000.0) + " sec\n");
                
//                SparseDoubleMatrix2D P = new SparseDoubleMatrix2D(transNormalized);
//                
//                System.out.print("zSum: " + P.viewColumn(0).zSum() + "\n");
//                P.viewColumn(0).assign(cern.jet.math.Functions.div(14)); 
//
//                System.out.println(P);
               
                m.add(i, RandomWalkWithRestart(P,proteinList.size(),i,a));
                
                //UNDO, get orig from transition an normalize
                //RestoreTransitionMatrix(i,transNormalized);                               
                for(int b=0;b<toBias.size();b++)
                {   
                    //double sum=0;
                    int currentNode = toBias.get(b);
                    P.viewColumn(currentNode).assign(Pt.viewColumn(currentNode));
                    P.viewColumn(currentNode).assign(cern.jet.math.Functions.div(Pt.viewColumn(currentNode).zSum()));

                }
                //System.out.println(P);
            }
            else
            {
                //SparseDoubleMatrix2D P = new SparseDoubleMatrix2D(transNormalized);
                //System.out.println(P);
                m.add(i, RandomWalkWithRestart(P,proteinList.size(),i,a));
            }
            //count++;
            //endTime = System.nanoTime();
            //System.out.println(count + " : " + ((endTime-startTime)/1000000000.0) + " sec\n");
                        
        }
        endTime = System.nanoTime();
        System.out.println("Elaspsed Time: " + ((endTime-startTime)/1000000000.0) + " sec\n");
        //Define set of clusters D ArrayList of an ArrayList of ints 
        //steps RRW 5-6
        ArrayList<ArrayList<Integer>> D = new ArrayList<ArrayList<Integer>>();
        System.out.println("Clustering with RepeatedRandomWalks..."); 
        for(int i=0;i<proteinList.size();i++)
        {
            //define Arraylist of clusters - starts empty
            ArrayList<ArrayList<Integer>> W = new ArrayList<ArrayList<Integer>>();
            //
            ArrayList<Integer> N = new ArrayList<Integer>();
            //add the current node to the first cluster
            N.clear();
            N.add(i);
            double previousWeight=0;
            //Create new cluster
            ArrayList<Integer> F = new ArrayList<Integer>();
            W.add(N);
            while(!W.isEmpty())
            { 
                ArrayList<Integer> C = new ArrayList<Integer>();
                C = W.remove(0); //Take out the current cluster to expand it
                DenseDoubleMatrix1D B = ClusterWalks(C, m);
                //get current (cluster size+1)th largest max node value and index
                double [] b = B.toArray();
                int c=0;
                double Bc=0;
                for(int j=0;j<b.length;j++)
                {
                    //if(i==vif)
                        //System.out.print("Is " + b[j] + " > " + Bc+ " and " + !C.contains(j) + '\n');
                    if((b[j] > Bc) && !C.contains(j) && (transition[j][C.get(0)] > 0)) //i!=j because we don't want the actual node
                    {
                        //if(i==vif)
                        //    System.out.print("--Yes " + b[j] + " > " + Bc+ " and " + !C.contains(j) + '\n');
                        c = j;
                        Bc = b[j];
                    }
                }
                //if(i==vif)
                //   System.out.print("Next largest " + Bc + '\n');
                //if within cutoff, add next largest value of B to the cluster
                if((Bc >= lamda*previousWeight))
                {
                    //if(i==vif)
                    //    System.out.print("ADDING protein " + proteinList.get(c) + " with score " + Bc + " >= (" + lamda + " * " +previousWeight + ") = " + lamda*previousWeight + '\n' + '\n');
                    F = C;
                    F.add(c);
                    W.add(F); //
                    previousWeight = Bc; //save for adding next node
                }
                else{
                    //if(i==vif)
                    //    System.out.print("NOT ADDING protein " + proteinList.get(c) + " with score " + Bc + " < (" + lamda + " * " +previousWeight + ") = " + lamda*previousWeight + '\n');
                }
                //if more than max cluster size, break
                if(!(F.size() < k))
                {
                    break;
                }
            }
/********************************************************************************/
            //CHECK MINIMUM CONNECTIONS
            //CHECK FOR DUPLICATE CLUSTERS AND MINIMUM SIZE
            //System.out.println("Checking for minimum size and duplicates...");
            boolean duplicate = false;
            if(F.size() > 1) //minimum size 3
            {
                //check to make sure not a total dup
                for(int j=0;j<D.size();j++)
                {
                    ArrayList cluster = new ArrayList();
                    cluster = D.get(j);
                    if(cluster.containsAll(F))
                    {
                        duplicate=true;
                        break;
                    }
                }
                //count the edges
                int connections=0;
                for(int j=0;j<F.size()-1;j++)
                {
                    for(int l=j;l<F.size()-1;l++)
                    {
                        if(transition[F.get(j)][F.get(l+1)] != 0)
                        {
                            connections++;
                        }
                    }
                }
                //System.out.println("Cluster of size "+F.size() + " and connections "+ connections);
                //if(!duplicate && (connections > F.size())) //USE LCC?
                if(!duplicate && (connections > 1)) //USE LCC?
                //if(!duplicate) //USE LCC?
                {
                    D.add(F); //Add F to the set of clusters
                }
            }
        }
        endTime = System.nanoTime(); 
        System.out.println("Elaspsed Time: " + ((endTime-startTime)/1000000000.0) + " sec\n");   
        
        System.out.println("Printing Clusters...");
        PrintClusters(D,proteinList);
        endTime = System.nanoTime(); 
        System.out.println("Elaspsed Time: " + ((endTime-startTime)/1000000000.0) + " sec\n");  
    }
    //Print Clusters to Console
       public static void PrintClusters(ArrayList<ArrayList<Integer>> D, ArrayList<String> proteinList){
        for(int i=0;i<D.size();i++)
        {
            int count = 0;
            for(int p =0;p<pathogenList.size();p++)
                if(D.get(i).contains(pathogenList.get(p)))
                {
                    count++;
                }
            
            System.out.print("\n");
            ArrayList<Integer> cluster = new ArrayList<Integer>();
            cluster = D.get(i);
            //double edges[] cluster.;
            ArrayList<Double> edges = new ArrayList<Double>();
            for(int j=0;j<cluster.size()-1;j++)
            {
                for(int l=j;l<cluster.size()-1;l++)
                {
                    if(transition[cluster.get(j)][cluster.get(l+1)] != 0)
                    {
                        edges.add(transition[cluster.get(j)][cluster.get(l+1)]);
                    }
                }
            }
            double out[] = new double[edges.size()];
            for(int d=0;d<edges.size();d++){
                out[d] = edges.get(d);
            }
            System.out.print(count + " " + cluster.size()+ " " + edges.size() + " "+ sum(out) + " " + average(out)+ " " + variance(out)+ " " + min(out)+ " " + max(out)+ " " + median(out)+ " ");
            for(int j=0;j<cluster.size();j++)
            {
                System.out.print(proteinList.get(Integer.parseInt(cluster.get(j).toString()))+" ");
            }
            
        }
        System.out.print("\n"); 
        System.out.print("\n");         
    }
    public static double sum(double a[]) {
        double sum = 0.0;
        for (double num : a)
            sum += num;
        return sum;
    }
    public static double average(double a[]) {
        double sum = 0.0;
        for (double num : a)
            sum += num;
        return sum/(a.length);
    }
    public static double variance(double a[]) {
        double sum = 0.0;
        for (double num : a)
            sum += num;
        double v = 0;
        for(double confidence : a)
            v += (sum/(a.length)-confidence)*(sum/(a.length)-confidence);
        return v/(a.length);
    }
    public static double min(double a[]) {
        double min = a[0];
        for(int i=0; i<a.length;i++)
        {
        if(a[i]<min)
            min = a[i];
        }
        return min;
    }
    public static double max(double a[]) {
        double max = a[0];
        for(int i=0; i<a.length;i++)
        {
        if(a[i]>max)
            max = a[i];
        }
        return max;
    }
    public static double mode(double a[]) {
        double maxValue=0.0;
        int maxCount=0;
        for (int i = 0; i < a.length; ++i) {
            int count = 0;
            for (int j = 0; j < a.length; ++j) {
                if (a[j] == a[i])
                    ++count;
            }
            if (count > maxCount) {
                maxCount = count;
                maxValue = a[i];
            }
        }
        return maxValue;
    }
    public static double median(double a[]) {
        Arrays.sort(a);  //Sort the array
        int mid = a.length / 2;
        if (a.length % 2 == 1) {
            return a[mid];
        } else {
            return (a[mid - 1] + a[mid])/ 2.0;
        }
    }
    //Prints the degrees of the vertices
    public static void StoreDegreeDistribution(double transition[][]){
       //System.out.println("Degree Distribution:");
        
       degree.clear();
       
       for(int i=0;i<transition.length;i++)
       {
           int count=0;
           for(int j=0;j<transition[i].length;j++)
           {
               if(transition[i][j] !=0)
                count++;
           }
           
           degree.add(count);
       }
    }  
    
    //Prints the degrees of the vertices
    public static void CalcKcore(ArrayList<Integer> d){

        while(Collections.min(degree) != 999999){
            
            ArrayList<Integer> kcore = new ArrayList<Integer>();
            int minD = Collections.min(degree);
            int node = degree.indexOf(minD);            
        
            while((Collections.min(degree) <= minD)){

                //add v to current core
                kcore.add(node);
                //cal degrees again (all or just subtract)
                for(int j=0;j<transition.length;j++)
                {
                    if((transition[node][j] !=0) && degree.get(j) != 999999)
                        degree.set(j, degree.get(j) - 1);
                }
                //remove v and all its edges            
                //degree.remove(node);
                degree.set(node, 999999);
                //minD = Collections.min(degree);
                node = degree.indexOf(Collections.min(degree));  
            }
            kCores.add(kcore);
            double edges=0.0;
            for(int i=0;i<kcore.size();i++)
            {   
                nodeCoreSize.set(kcore.get(i), kcore.size());
                
                for(int j=i;j<kcore.size()-1;j++)
                {
                    if((transition[kcore.get(i)][kcore.get(j+1)] !=0))
                        edges++;
                }

            }  
            for(int i=0;i<kcore.size();i++)
            {   
                nodeCoreDensity.set(kcore.get(i),edges/((kcore.size()*(kcore.size()-1))/2.0));
            }              
                        
            
            
        }
    }     


}