package wmoFSS;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import wmoFSS.MOFish;
import General.Run;
import java.util.List;

public class MOSchool {

	public double[][][] history;
	private ArrayList<MOFish> fishes;
	private double stepInd;
	private double stepVol;
	private double alpha;
	
	public double[] fMin;
	public double[] fMax;
	
	private HashMap <Integer, ArrayList<MOFish>> clusters;
	private ArrayList<MOFish> finalPareto;
	public moFitnessFunction fitnessFunction;
	private double IGD;
	
	private ArrayList<MOFish> extremePoints;
	private List<Double> intercepts;

	
	public MOSchool() throws FileNotFoundException{
		
		fitnessFunction = new moFitnessFunction();
				
		this.fishes = new ArrayList<MOFish>();
		this.history = new double[Run.NUMBER_OF_ITERATIONS][Run.NUMBER_OF_FISH][Run.NUMBER_OF_OBJECTIVES];
		
		finalPareto = new ArrayList<MOFish>();
		
		fMin = new double[Run.NUMBER_OF_OBJECTIVES];
		fMax = new double[Run.NUMBER_OF_OBJECTIVES];
		
		
		clusters = new HashMap<Integer, ArrayList<MOFish>>((fitnessFunction.NUMBER_OF_REFPOINTS));
		
		for(int i=0; i<fitnessFunction.NUMBER_OF_REFPOINTS; i++){
			ArrayList<MOFish> fishInCluster = new ArrayList<MOFish>();
			clusters.put(i, fishInCluster);
		}
				
		for(int i = 0; i<fMin.length; i++) {
			fMin[i] = Double.POSITIVE_INFINITY;
			fMax[i] = Double.NEGATIVE_INFINITY;
		}
		
		for(int i = 0; i<Run.NUMBER_OF_FISH; i++){
			fishes.add(new MOFish(fitnessFunction));			
		}
				
		stepInd = Run.INITIAL_IND_STEP;
		stepVol = Run.INITIAL_VOL_STEP;
		
		extremePoints = new ArrayList<MOFish>();
		intercepts = new ArrayList<Double>();
		
		for(int i=0; i<fMin.length; i++) {
			this.fMin[i] = 0.0;			
		}
		
	}
	
	public void runSearch(){
		//Initial Clustering
		updateBoundaries();
		feeding();					// Clusterização apenas no começo da busca
		
		if (Run.GEOMETRIC_CLUSTERING){
			clusteringGeometric();  // Clusterização apenas no começo da busca
		} else {
			clustering();			// Clusterização apenas no começo da busca
		}
		
		splitInClusters();
		
		if(Run.INITIAL_CLUSTERING){ //Only Initial Clustering
			
			for(int i=0; i<Run.NUMBER_OF_ITERATIONS;i++){
				updateBoundaries();
				feeding();
				individualMovement(i);
				updateBoundaries();
				feeding();
				leadersDefine();
				collectiveInstintiveMovement();
				collectiveVolitiveMovement();
				updateSteps(i);
				//updateHistory(i);
				
				updateFitnesses();
				
				//paretoSorting();
				
				updateHistory(i);
				
			}
		}else{ //Clustering Every Iteration
			for(int i=0; i<Run.NUMBER_OF_ITERATIONS;i++){
				updateBoundaries();
				feeding();
				individualMovement(i);
				updateBoundaries();
				//feeding();
				if (Run.GEOMETRIC_CLUSTERING){
					clusteringGeometric();  // Clusterização apenas no começo da busca
				} else {
					clustering();			// Clusterização apenas no começo da busca
				}
				
				feeding();		// Para calcular o novo peso agregado, já que este depende do cluster em que o peixe está
				splitInClusters();
				
				collectiveInstintiveMovement();
				collectiveVolitiveMovement();
				updateSteps(i);
				//updateHistory(i);
				
				updateFitnesses();
				
				//paretoSorting();
				
				updateHistory(i);
				
			}
		}
		
		updateBoundaries();
		feeding();
		leadersDefine();
		//returnLeaders();
		paretoSortingClusters();
		//System.out.println("Final Pareto - Positions");
		//printListOfFishesPositions(finalPareto);
		System.out.println("Final Pareto - Fitness");
		printListOfFishesFitness(finalPareto);
		
		//System.out.println("Final Pareto - Peso");
		//printListOfFishesWeights(finalPareto);
		
		IGD = calculateIGD(fitnessFunction.paretoOptimalFront, finalPareto);
		System.out.println(IGD);
		//double HV = calculateHV(fitnessFunction.getReferencePointHV(), finalPareto);
		//System.out.println(HV);
	}
	
	private void individualMovement(int iteration) {
		
		for(MOFish e:fishes){
			
			e.moveFishIndividual(stepInd, iteration, fMin, fMax, intercepts, alpha);
						
		}
	}
	
	private void updateBoundaries() {		// Atualiza os máximos e mínimos dos valores de fitness encontrados
		
		if (Run.ASF_NORMALIZATION) {
			
			extremePoints = findExtremePoints(fishes);
			intercepts =  constructHyperplane(fishes, this.extremePoints);
		
		}else{
		
			double[] maxFToUpdate = getMaxF();
			//double[] minFToUpdate = getMinF();
			
			fMax = maxFToUpdate;
						
		}
		
	}
	
	private double[] getMinF() {
		
		double[] minF = new double[Run.NUMBER_OF_OBJECTIVES];
		double[] fishFitness = new double[Run.NUMBER_OF_OBJECTIVES];
		
		for(int i=0; i<minF.length; i++) {
			minF[i] = Double.POSITIVE_INFINITY;
		}
		
		for(MOFish e:fishes){
			fishFitness = e.getFitness();
			for (int i=0; i<minF.length; i++) {
				if( Double.compare(fishFitness[i], minF[i]) < 0) {
					minF[i] = fishFitness[i];
				}
			}
		}
		
		return minF;
	}

	private double[] getMaxF() {
		
		double[] maxF = new double[Run.NUMBER_OF_OBJECTIVES];
		double[] fishFitness = new double[Run.NUMBER_OF_OBJECTIVES];
		
		for(int i=0; i<maxF.length; i++) {
			maxF[i] = Double.NEGATIVE_INFINITY;
		}
		
		for(MOFish e:fishes){
			fishFitness = e.getFitness();
			for (int i=0; i<maxF.length; i++) {
				if(Double.compare(fishFitness[i], maxF[i]) > 0){
					maxF[i] = fishFitness[i];
				}
			}
		}
		
		return maxF;
	}
	
	private void feeding(){
		
		//double[] lastWeights = new double[Run.NUMBER_OF_FISH];
		
		if(Run.NORMALIZATIONONOFF){
			
			if (Run.ASF_NORMALIZATION) {
				
				for(MOFish e:fishes){
					e.feedFish2(this.intercepts, fMin);
				}
				
			}else{
				
				for(MOFish e:fishes){
					e.feedFish(fMax, fMin);
				}
				
			}
			
		}else{
			
			for(MOFish e:fishes){
				e.feedFishNotNormalizing();
			}
			
		}
			
	}
	
	private void clustering() {		// Separa os peixes em clusters de acordo com a distância euclidiana mínima para os pontos de referência
		
		// Cleaning clusters hash map
		Set<Integer> keysSet = clusters.keySet();
		
		for(Integer e:keysSet){
			ArrayList<MOFish> fishInCluster = clusters.get(e);
			fishInCluster.clear();
		}
		
		// New clustering
			for(MOFish e:fishes) {		// Neste laço cada a distância entre cada peixe e cada ponto de referência é calculada
										// O peixe é colocado no cluster do ponto de referência mais próximo a ele
				int closerRefPoint = 0;
				double smallerDistance = Double.POSITIVE_INFINITY;
				
				for(int i=0; i<fitnessFunction.getReferencePoints().size(); i++){
					double distance = perpendicularDistance(fitnessFunction.getReferencePoints().get(i), e.getWeight());
					if (Double.compare(distance, smallerDistance) < 0) {
						smallerDistance = distance;
						closerRefPoint = i;
					}
				}
				
				e.setRefPointIndex(closerRefPoint);		// Atualiza o ponto de referência ao qual o peixe "e" está associado
				
				ArrayList<MOFish> fishInCluster = clusters.get(closerRefPoint);  // Coloca o peixe "e" no cluster de "closerRefPoint"
				fishInCluster.add(e);
				//clusters.remove(closerRefPoint);
				//clusters.put(closerRefPoint, fishInCluster);
				
			}
	}
	
	private void clusteringGeometric() {
		     // Separa os peixes em clusters de acordo com a distância euclidiana mínima para os pontos de referência
            
            double fishesPerCluster = Math.round(Run.NUMBER_OF_FISH/fitnessFunction.NUMBER_OF_REFPOINTS);
            
            ArrayList<MOFish> fishesForClustering = new ArrayList<MOFish>();

            for(MOFish a : fishes){
                   fishesForClustering.add(a);
            }            
           
            HashMap<Double, MOFish> fishesList;
           
            ArrayList<double[]> refPointsList = fitnessFunction.getReferencePoints();
           
            // Cleaning clusters hash map
            Set<Integer> keysSet = clusters.keySet();
           
            for(Integer e:keysSet){
                   ArrayList<MOFish> fishInCluster = clusters.get(e);
                   fishInCluster.clear();
            }
           
            // New clustering
                   for(int j=0; j<refPointsList.size(); j++) {
                        
                	   for(int i=0; i<fishesPerCluster; i++){
                               
                		   double smallerDistance = Double.POSITIVE_INFINITY;
                           fishesList = new HashMap<Double, MOFish>();
                           int indexToRemove = -1;
                           
                           for(int k=0; k<fishesForClustering.size(); k++){
                        	   double distance = perpendicularDistance(refPointsList.get(j), fishesForClustering.get(k).getWeight());
                               
                        	   fishesList.put(distance, fishesForClustering.get(k));
                               
                        	   if (Double.compare(distance, smallerDistance) < 0) {
                            	   smallerDistance = distance;
                            	   indexToRemove = k;
                               }
                           }
                               
                           MOFish FishToSet = fishesList.get(smallerDistance);
                           fishesForClustering.remove(indexToRemove);
                           
                           FishToSet.setRefPointIndex(j);
                           
                        	                               
                           ArrayList<MOFish> fishInCluster = clusters.get(j);
                           fishInCluster.add(FishToSet);
                               
                	   }
                	   
                   }
                   
                   if(fishesForClustering.isEmpty() == false) {
                   	
                	   for(MOFish e:fishesForClustering) {
                		   
                		   int closerRefPoint = 0;
                		   double smallerDistance = Double.POSITIVE_INFINITY;
                		   
                		   for(int i=0; i<fitnessFunction.getReferencePoints().size(); i++){
                			   double distance = perpendicularDistance(fitnessFunction.getReferencePoints().get(i), e.getFitness());
                			   if (Double.compare(distance, smallerDistance) < 0) {
                				   smallerDistance = distance;
                				   closerRefPoint = i;                				   
                			   }
                			   
                		   }
                		   
                		  e.setRefPointIndex(closerRefPoint);
                		  //System.out.println("Ref Point Index " + closerRefPoint);
                		  //printArray(e.getFitness());
                		  ArrayList<MOFish> fishInCluster = clusters.get(closerRefPoint);
                		  fishInCluster.add(e);
                		  //clusters.remove(closerRefPoint);
                		  //clusters.put(closerRefPoint, fishInCluster);
                		  
                	   }
                	   
                   }
                   
	}
	
	private double perpendicularDistance(double[] refPoint, double[] solution) {
		double d1 = 0.0;
		double d2 = 0.0;
		
		// Cálculo de d1
		double num1 = 0.0;
		double den1 = 0.0;
		for (int i=0; i<Run.NUMBER_OF_OBJECTIVES; i++){
			num1 += solution[i]*refPoint[i];		// Produto interno entre os vetores
			den1 += Math.pow(refPoint[i], 2);			// Cálculo da norma de refPoint
		}
		den1 = Math.sqrt(den1);							// Cálculo da norma de refPoint
		d1 = Math.abs(num1)/den1;
		
		// Cálculo de d2
		for (int i=0; i<Run.NUMBER_OF_OBJECTIVES; i++){
			double aux1 = (d1*refPoint[i])/den1;			// Cálculo de w - (d1*lambda)/(norma de lambda)
			double aux2 = solution[i] - aux1;				// Cálculo de w - (d1*lambda)/(norma de lambda)
			d2 += Math.pow(aux2, 2);
		}
		d2 = Math.sqrt(d2);
		
		return d2;
	}
	
	private void splitInClusters() {	// Após os processo de clusterização, cada  peixe deve receber a ArrayList que corresponde a seu cluster 				
		for(MOFish e:fishes){
			e.setFishCluster(clusters.get(e.getRefPointIndex()));
		}
		
		leadersDefine();
	}
	
	private void leadersDefine(){
		
		double minWeight;
		Set<Integer> keysSet = clusters.keySet();
		ArrayList<MOFish> clusterToDefineLeader;
		MOFish leaderToSet = new MOFish(fitnessFunction);
		
		for(Integer e:keysSet){
			
			clusterToDefineLeader = clusters.get(e);
			minWeight = Double.POSITIVE_INFINITY;
			
			for(MOFish a:clusterToDefineLeader){
				
				a.resetLeader();
				if (Double.compare(a.getAggregatedWeight(), minWeight) < 0) {
					minWeight = a.getAggregatedWeight();
					leaderToSet = a;
				}
				
			}
			
			leaderToSet.setLeader();
			
			
		}
		
	}
	
	
	private void collectiveInstintiveMovement() {
		
		for(MOFish e:fishes){
			e.calculateI();
		}
		
		for(MOFish a:fishes){
			a.moveFishColInstinctive();
		}
		
	}
	
	private void collectiveVolitiveMovement() {
		
		for(MOFish e:fishes){
			e.calculateB();
		}
		
		for(MOFish a:fishes){
			a.moveFishVolitive(stepVol);
		}		
	}
	
	private void updateSteps(int currentIteration){
		stepInd = stepInd - (Run.INITIAL_IND_STEP - Run.FINAL_IND_STEP) / Run.NUMBER_OF_ITERATIONS;
		stepVol = stepVol - (Run.INITIAL_VOL_STEP - Run.FINAL_VOL_STEP) / Run.NUMBER_OF_ITERATIONS;
		//stepInd = Run.INITIAL_IND_STEP*Math.exp(-0.001*currentIteration);
		//stepVol = Run.INITIAL_VOL_STEP*Math.exp(-0.007*currentIteration);
		
		alpha = alpha - (Run.INITIAL_ALPHA)/Run.NUMBER_OF_ITERATIONS;
		
	}
		
	
	private double euclidianDistance(double[] refPoint, double[] fishFit){
		double distance = 0;
		double diff = 0;
		for (int i=0; i<fishFit.length; i++) {
			diff = refPoint[i] - fishFit[i];
			distance = distance + Math.pow(diff, 2); 
		}
		
		distance = Math.sqrt(distance);
		
		return distance;
	}
	
	
	// CALCULA A MENOR DISTÂNCIA ENTRE UMA SOLUÇÃO DO PARETO OPTIMAL FRONT E AS SOLUÇÕES DO PARETO ENCONTRADO
	private double minimumDistance(double[] paretoOptimalFrontMember, ArrayList<MOFish> finalPareto) {
		double minimum = Double.POSITIVE_INFINITY;

		for (int i = 0; i < finalPareto.size(); i++) {
			minimum = Math.min(minimum, euclidianDistance(paretoOptimalFrontMember, finalPareto.get(i).getFitness()));
		}
		
		return minimum;
	}
	
	
	// CÁLCULO DO IGD
	private double calculateIGD(ArrayList<double[]> paretoOptimalFront, ArrayList<MOFish> finalPareto){
		double IGD = 0;
		
		for (int i = 0; i < paretoOptimalFront.size(); i++) {
			IGD += Math.pow(minimumDistance(paretoOptimalFront.get(i), finalPareto), 2);
		}
		
		IGD = Math.sqrt(IGD)/paretoOptimalFront.size();
		
		return IGD;
	}
	
	private double minimumDistance2(double[] paretoOptimalFrontMember, ArrayList<double[]> finalPareto) {
		double minimum = Double.POSITIVE_INFINITY;

		for (int i = 0; i < finalPareto.size(); i++) {
			minimum = Math.min(minimum, euclidianDistance(paretoOptimalFrontMember, finalPareto.get(i)));
		}
		
		return minimum;
	}
	
	
	private double calculateIGD2(ArrayList<double[]> paretoOptimalFront, ArrayList<double[]> finalPareto){
		double IGD = 0;
		
		for (int i = 0; i < paretoOptimalFront.size(); i++) {
			IGD += Math.pow(minimumDistance2(paretoOptimalFront.get(i), finalPareto), 2);
		}
		
		IGD = Math.pow(IGD, 1/2)/paretoOptimalFront.size();
		
		return IGD;
	}
	
	public double getIGD(){
		return IGD;
	}
	

	private void updateFitnesses() {
		
		for(MOFish e:fishes){
			double[] newFitness = new double[Run.NUMBER_OF_OBJECTIVES];
			newFitness = e.calculateFitness(e.getPosition());
			e.updateFitness(newFitness);
		}
	}
	
	// private void updateWeights() {
		
	//	for(MOFish e:fishes){
	//		e.updateWeight(e.feedFish(fMax, fMin));
	//	}
	// }
	
	private void paretoSorting() {
		finalPareto.clear();
		for (MOFish e:fishes) {
			boolean nonDominated = true;
			for (MOFish f:fishes) {
				int dominate = e.paretoCompare(f.getFitness());
				if (dominate == 2) {		// Caso qualquer solução domine "e", nonDominated se torna falso  
					nonDominated = false;
				}
			}
			if (nonDominated == true) {
				finalPareto.add(e);
			}
		}
	}
	
	private void paretoSortingClusters() {
		finalPareto.clear();
		Set<Integer> keysSet = clusters.keySet();
		ArrayList<MOFish> clusterToSort;
		
		
		for(Integer e:keysSet){
			
			clusterToSort = clusters.get(e);
			
			for(MOFish a:clusterToSort){
				boolean nonDominated = true;
				for (MOFish f:fishes) {
					int dominate = a.paretoCompare(f.getFitness());
					if (dominate == 2) {		// Caso qualquer solução domine "e", nonDominated se torna falso  
						nonDominated = false;
					}
				}
				
				if (nonDominated == true) {
					finalPareto.add(a);
				}

			}
				
		}			
	}
	
	private void returnLeaders() {
		finalPareto.clear();
		Set<Integer> keysSet = clusters.keySet();
		ArrayList<MOFish> clusterToGetLeader;
		
		
		for(Integer e:keysSet){
			
			clusterToGetLeader = clusters.get(e);
			
			for(MOFish a:clusterToGetLeader){
				if (a.getIsLeader()) {
					finalPareto.add(a);
				}	
			}
		}
	}	

	
	private double ASF(MOFish fish, int index) {
		
		double max_ratio = Double.NEGATIVE_INFINITY;
		
		for (int i = 0; i < Run.NUMBER_OF_OBJECTIVES; i++) {
			double weight = (index == i) ? 1.0 : 0.000001;
			max_ratio = Math.max(max_ratio, fish.getFitness()[i]/weight);
		}
		
		return max_ratio;
	}
	
	private ArrayList<MOFish> findExtremePoints(ArrayList<MOFish> fishes) {
		ArrayList<MOFish> extremePoints = new ArrayList<>();
		MOFish min_indv = null;
		for (int f=0; f < Run.NUMBER_OF_OBJECTIVES; f+=1)
		{
			double min_ASF = Double.MAX_VALUE;	
			for (MOFish s : fishes) { // 
				double asf = ASF(s, f);
				if (Double.compare(asf, min_ASF) < 0) {
					min_ASF = asf;
					min_indv = s;
				}
			}
			
			extremePoints.add(min_indv);
		}
		return extremePoints;
	}
	
	public List<Double> gaussianElimination(List<List<Double>> A, List<Double> b) {
		List<Double> x = new ArrayList<>();

	    int N = A.size();
	    for (int i=0; i<N; i+=1)
	    {
	    	A.get(i).add(b.get(i));
	    }

	    for (int base=0; base<N-1; base+=1)
	    {
	        for (int target = base+1; target<N; target+=1)
	        {
	            double ratio = A.get(target).get(base)/A.get(base).get(base);
	            for (int term=0; term < A.get(base).size(); term += 1)
	            {
	                A.get(target).set(term, A.get(target).get(term) - A.get(base).get(term)*ratio);
	            }
	        }
	    }

	    for (int i = 0; i < N; i++)
	    	x.add(0.0);
	    
	    for (int i=N-1; i>=0; i-=1)
	    {
	        for (int known=i+1; known<N; known+=1)
	        {
	            A.get(i).set(N, A.get(i).get(N) - A.get(i).get(known)*x.get(known));
	        }
	        x.set(i, A.get(i).get(N)/A.get(i).get(i));
	    }
		return x;
	}
	
	public List<Double> constructHyperplane(ArrayList<MOFish> fishes, ArrayList<MOFish> extreme_points) {
		// Check whether there are duplicate extreme points.
		// This might happen but the original paper does not mention how to deal with it.
		boolean duplicate = false;
		for (int i=0; !duplicate && i< extreme_points.size(); i+=1)
		{
			for (int j=i+1; !duplicate && j<extreme_points.size(); j+=1)
			{
				duplicate = extreme_points.get(i).equals(extreme_points.get(j));
			}
		}

		List<Double> intercepts = new ArrayList<>();
		
		if (duplicate) // cannot construct the unique hyperplane (this is a casual method to deal with the condition)
		{
			for (int f=0; f<Run.NUMBER_OF_OBJECTIVES; f+=1)
			{
				// extreme_points[f] stands for the individual with the largest value of objective f
				intercepts.add(extreme_points.get(f).getFitness()[f]);
			}
		}
		else
		{
			// Find the equation of the hyperplane
			List<Double> b = new ArrayList<>(); //(pop[0].objs().size(), 1.0);
			for (int i =0; i < Run.NUMBER_OF_OBJECTIVES;i++)
				b.add(1.0);
			
			List<List<Double>> A=new ArrayList<>();
			for (MOFish s : extreme_points)
			{
				List<Double> aux = new ArrayList<>();
				for (int i = 0; i < Run.NUMBER_OF_OBJECTIVES; i++)
					aux.add(s.getFitness()[i]);
				A.add(aux);
			}
			List<Double> x = gaussianElimination(A, b);
		
			// Find intercepts
			for (int f=0; f<Run.NUMBER_OF_OBJECTIVES; f+=1)
			{
				intercepts.add(1.0/x.get(f));
				
			}
		}
		return intercepts;
	}
	
	private double calculateHV(double[] referencePoint, ArrayList<MOFish> finalPareto) {

		double[] randomResponseToCheck = new double[Run.NUMBER_OF_OBJECTIVES];

		int inVolumeCounter = 0;		
		double ans;
		int test;
		boolean check;

		for(int i=0; i<Run.numberOfRandomResponses; i++){

			randomResponseToCheck = generateRandomResponse(referencePoint);
			//printArray(randomResponseToCheck);

			check = true;

			for(MOFish e:finalPareto){
				
				test = e.paretoCompare(randomResponseToCheck);
				//System.out.println(test);
				
				if(test == 2) {
					check = false;
				}

			}

			if (check) {

				inVolumeCounter++;

			}

		}
		
		ans = inVolumeCounter/Run.numberOfRandomResponses;

		return ans;

	}
	
	private double[] generateRandomResponse(double[] referencePoint){

		double[] randomResponseToReturn = new double[referencePoint.length];

		double randomNumber;

		for(int i=0; i<referencePoint.length; i++){

			randomNumber = Run.generator.nextDouble();

			randomResponseToReturn[i] = referencePoint[i]*randomNumber;

		}

		return randomResponseToReturn;

	}

	private void updateHistory(int currentIteration){
		
		int i = 0;
		
		for (MOFish e: fishes){
			this.history[currentIteration][i] = e.getFitness();
			i++;
		}
	}
	
	public ArrayList<double[]> getFinalPareto() {
		ArrayList<double[]> finalParetoFitness = new ArrayList<double[]>();
		
		for (MOFish fish:finalPareto) {
			finalParetoFitness.add(fish.getFitness());
		}
		
		return finalParetoFitness;
	}
	
	public static void printListOfFishesWeights(ArrayList<MOFish> listToPrint){
		//System.out.println("Lista = " + listToPrint);
		for(MOFish fish: listToPrint){
			for(int i=0; i< fish.getWeight().length; i++) {
			System.out.println("Component[" + i + "] = " + fish.getWeight()[i]);
			}
		}
	}

	public static void printListOfFishesFitness(ArrayList<MOFish> listToPrint){
		//System.out.println("Lista = " + listToPrint);
		for(MOFish fish: listToPrint){
			for(int i=0; i< fish.getFitness().length; i++) {
			System.out.println("Component[" + i + "] = " + fish.getFitness()[i]);
			}
		}
	}
	
	public static void printListOfFishesPositions(ArrayList<MOFish> listToPrint) {
		for(MOFish fish: listToPrint){
			for(int i=0; i< fish.getPosition().length; i++) {
			System.out.println("Component[" + i + "] = " + fish.getPosition()[i]);
			}
		}
	}
	
	public static void printArray(double[] arrayToPrint){
		System.out.println("Array = " + arrayToPrint);
			for(int i=0; i< arrayToPrint.length; i++) {
			System.out.println("Component[" + i + "] = " + arrayToPrint[i]);
			}
	}

	public static void printHashMap(HashMap <Integer, ArrayList<MOFish>> hashToPrint){
		System.out.println("HashMap = " + hashToPrint);
			for(int i=0; i< hashToPrint.size(); i++) {
				printListOfFishesFitness(hashToPrint.get(i));
			}
	}
	
}