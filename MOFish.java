package wmoFSS;

import java.util.ArrayList;
import java.util.List;

//import org.moeaframework.core.Solution;

import General.Run;

public class MOFish {
	
	private double[] position;
	private double[] lastPosition;
	private double[] fitness;
	private double[] lastFitness;
	private double[] lastWeight;
	private double[] weight;
	private double aggregatedWeight;
	private double lastAggregatedWeight;
	private double[] I;
	private double[] B;
	private double sense;
	private boolean leader;
	private int refPointIndex;
	private ArrayList<MOFish> fishCluster;
	private moFitnessFunction moFitnessFunction;
		
	public MOFish(moFitnessFunction fitnessFunction){
		
		position = new double[Run.DIMENSIONS];
		lastPosition = new double[Run.DIMENSIONS];
		fitness = new double [Run.NUMBER_OF_OBJECTIVES];
		lastFitness = new double [Run.NUMBER_OF_OBJECTIVES];
		lastWeight = new double [Run.NUMBER_OF_OBJECTIVES];
		weight = new double [Run.NUMBER_OF_OBJECTIVES];
		I = new double[Run.DIMENSIONS];
		B = new double[Run.DIMENSIONS];
		fishCluster = new ArrayList<MOFish>();
		moFitnessFunction = fitnessFunction;
				
		position = generateRandomPosition();
		fitness = calculateFitness(position);
		
		updatePosition(position);
		updateFitness(fitness);
	}

	public void updateWeight(double[] newWeight) {
		
		lastWeight = this.weight;
		weight = newWeight;
		
	}

	private void updatePosition(double[] newPosition) {
				
		lastPosition = this.position;
		position = newPosition;
		
	}

	public void updateFitness(double[] newFitness) {
		lastFitness = this.fitness;
		fitness = newFitness;
	}

	public double[] calculateFitness(double[] position) {
		
		double[] newFitness = new double[Run.NUMBER_OF_OBJECTIVES];
		newFitness = moFitnessFunction.calculateFitness(position);
		return newFitness;
	}
	
	public double[] getFitness() {
		return this.fitness;
	}

	public double[] getLastFitness() {
		return this.lastFitness;
	}
	
	public double[] getWeight() {
		return this.weight;
	}
	

	public double[] getPosition() {
		return position.clone();
	}
	
	public double[] getLastPosition() {
		return lastPosition.clone();
	}

	public void setFishCluster(ArrayList<MOFish> fishCluster) {
		this.fishCluster = fishCluster;
		//printListOfFishes(this.fishCluster);
	}

	private double aggregate(double[] weightVector){
		
		double[] refPoint = moFitnessFunction.getReferencePoints().get(this.refPointIndex);
		
		double theta = Run.theta;
		double d1 = 0;
		double d2 = 0;
		
		// Cálculo de d1
		double num1 = 0;
		double den1 = 0;
		for (int i=0; i<Run.NUMBER_OF_OBJECTIVES; i++){
			num1 += weightVector[i]*refPoint[i];		// Produto interno entre os vetores
			den1 += Math.pow(refPoint[i], 2);			// Cálculo da norma de refPoint
		}
		den1 = Math.sqrt(den1);							// Cálculo da norma de refPoint
		d1 = Math.abs(num1)/den1;
		
		// Cálculo de d2
		for (int i=0; i<Run.NUMBER_OF_OBJECTIVES; i++){
			double aux1 = (d1*refPoint[i])/den1;			// Cálculo de w - (d1*lambda)/(norma de lambda)
			double aux2 = weightVector[i] - aux1;				// Cálculo de w - (d1*lambda)/(norma de lambda)
			d2 += Math.pow(aux2, 2);
		}
		d2 = Math.sqrt(d2);
		
		double PBI = d1 + theta*d2;
		return PBI;
	}
	
	public double getAggregatedWeight() {	//recebe o vetor de pesos e retorna valor agregado
		
		double aggregatedWeight = aggregate(this.weight);
		return aggregatedWeight;
	}

	public double getLastAggregatedWeight() {
		lastAggregatedWeight = aggregate(this.lastWeight);
		return lastAggregatedWeight;
	}
	

	private double[] generateRandomPosition() {
		
		double xMin = Run.LOWEST_SPACE_BOUNDARY;
		double xMax = Run.HIGHEST_SPACE_BOUNDARY;
		double[] generatedPosition = new double[Run.DIMENSIONS];
		double randNum;
				
		for(int i=0;i<generatedPosition.length;i++){
			randNum = Run.generator.nextDouble();
			generatedPosition[i] = randNum * (xMax - xMin) + xMin;
		}
		
		return generatedPosition;
	}		
	
	public void moveFish(double[] displacement){
		double[] actualPosition = position;
		double[] finalPosition = new double[Run.DIMENSIONS];
		
		for(int i=0;i<displacement.length;i++){
			finalPosition[i] = actualPosition[i] + displacement[i];
		}
		
		this.updatePosition(bounderingControl(finalPosition));
		
	}
	
	
	private double[] findNewPosition(double[] displacement){
		
		double[] actualPosition = position;
		double[] finalPosition = new double[Run.DIMENSIONS];
		
		for(int i=0; i<displacement.length; i++){
			finalPosition[i] = actualPosition[i] + displacement[i];
		}
		
		return (bounderingControl(finalPosition));
		
	}
	
	private double[] bounderingControl(double[] positionToCheck){
		
		double[] toCorrect = positionToCheck;
		double xMin = Run.LOWEST_SPACE_BOUNDARY;
		double xMax = Run.HIGHEST_SPACE_BOUNDARY;
		
		for(int i=0;i<toCorrect.length;i++){
			
			if(Double.compare(toCorrect[i], xMax) > 0){ 
				toCorrect[i] = xMax;
			}else if(Double.compare(toCorrect[i], xMin) < 0){
				toCorrect[i] = xMin;
			}			
		}
		
		return toCorrect;
	}
		
	public int paretoCompare(double[] newFitness){
		
		double[] currentFitness = this.fitness;
		boolean dominate1 = false;
		boolean dominate2 = false;

		for (int i = 0; i < Run.NUMBER_OF_OBJECTIVES; i++) {
			if (Double.compare(newFitness[i], currentFitness[i]) > 0) {
				dominate1 = true;
					
				if (dominate2) {
					return 0;
				}
			} else if (Double.compare(currentFitness[i], newFitness[i]) > 0) {
				dominate2 = true;

				if (dominate1) {
					return 0;
				}
			}
		}

		if (dominate1 && !dominate2) {
			return 1;
		} else if (!dominate1 && dominate2) {
			return 2;
		} else {
			return 0;
		}
	}

	public void moveFishIndividual(double step, int iteration, double[] fMin, double[] fMax, List<Double> intercepts, double alphaLIN) {		// Já alterado para minimização
		
		double randNum;
		double[] newFitness;
		double[] displacement = new double[Run.DIMENSIONS];
		double[] candidatePosition = new double[Run.DIMENSIONS];
		double newAggregatedWeight;
		
		//SAR
		
		double alpha;
		
		switch(Run.sarMode){
		
			case 0: //No Worsening Allowance
				alpha = -1.0;
				break;
			case 1: //Exponential Decay Mode
				alpha = 0.8*Math.exp(-0.007*iteration);
				break;
			case 2: //Linear Decay Mode
				alpha = 1 - iteration/Run.NUMBER_OF_ITERATIONS;
				break;
			case 3: //Total Worsening Allowance
				alpha = 2;
				break;
			default:
				alpha = -1.0;	
		}
		
		//double alpha = alphaLIN;
		double randSAR = Run.generator.nextDouble();
		
		for(int i=0;i<displacement.length;i++){
			randNum = Run.generator.nextDouble();
			displacement[i] = 2*randNum*step-step;
		}
		
		candidatePosition = findNewPosition(displacement);
		newFitness = moFitnessFunction.calculateFitness(candidatePosition);

		if(Run.NORMALIZATIONONOFF){
			
			if (Run.ASF_NORMALIZATION) {
				
				newAggregatedWeight=this.testFeedFish2(intercepts, fMin, newFitness);
				
			}else{
				
				newAggregatedWeight=this.testFeedFish(fMax, fMin, newFitness);
				
			}
			
		}else{
			
			newAggregatedWeight=this.testFeedFishNotNormalizing(newFitness);
			
		}		
		
		if (Double.compare(newAggregatedWeight, this.aggregatedWeight) < 0) {
			updatePosition(candidatePosition);
			updateFitness(newFitness);
		} else {
			randNum = Run.generator.nextDouble();
			if (Double.compare(alpha, randSAR) > 0) {
				updatePosition(candidatePosition);
				updateFitness(newFitness);
			}
		}
		
	}
	
	public void feedFish(double[] maxF, double[] minF) {
		
		this.lastWeight = this.weight;
		
		for (int i=0; i<Run.NUMBER_OF_OBJECTIVES; i++){
			weight[i] = (fitness[i] - minF[i])/(maxF[i] - minF[i]);
		}
		this.lastAggregatedWeight = this.aggregatedWeight;
		this.aggregatedWeight = this.aggregate(weight);

	}
	
	
	public void feedFish2(List<Double> intercepts, double[] ideal_point) {
		
		
		for (int f = 0; f < Run.NUMBER_OF_OBJECTIVES; f++) {
			if (Math.abs(intercepts.get(f)-ideal_point[f])> 10e-10)
			{
				weight[f] = (fitness[f] - ideal_point[f]) / (intercepts.get(f) - ideal_point[f]);
			}
			else
			{
				weight[f] = (fitness[f] - ideal_point[f]) / (10e-10);
			}
				
		}
		
		this.lastAggregatedWeight = this.aggregatedWeight;
		this.aggregatedWeight = this.aggregate(weight);

	}
	
	public double testFeedFish(double[] maxF, double[] minF, double[] fitnessToTest) {
		
		double[] testWeight = new double[Run.NUMBER_OF_OBJECTIVES];
		double newAggregatedWeight;
		
		for (int i=0; i<Run.NUMBER_OF_OBJECTIVES; i++){
			testWeight[i] = (fitnessToTest[i] - minF[i])/(maxF[i] - minF[i]);
		}
		
		newAggregatedWeight = this.aggregate(testWeight);
		
		return newAggregatedWeight;
	}
	
	
	public double testFeedFish2(List<Double> intercepts, double[] ideal_point, double[] fitnessToTest) {
		
		double[] testWeight = new double[Run.NUMBER_OF_OBJECTIVES];
		double newAggregatedWeight;
		
		for (int f = 0; f < Run.NUMBER_OF_OBJECTIVES; f++) {
			if (Math.abs(intercepts.get(f) - ideal_point[f])> 10e-10)
			{
				testWeight[f] = (fitnessToTest[f] - ideal_point[f]) / (intercepts.get(f) - ideal_point[f]);
			}
			else
			{
				testWeight[f] = (fitnessToTest[f] - ideal_point[f]) / (10e-10);
			}
				
		}
		
		newAggregatedWeight = this.aggregate(testWeight);
		
		return newAggregatedWeight;

	}
	
	public double testFeedFishNotNormalizing(double[] fitnessToTest) {
		
		double newAggregatedWeight;
		
		newAggregatedWeight = this.aggregate(fitnessToTest);
		
		return newAggregatedWeight;
	}
	

	public void calculateI(){    // Recebe um peixe, seu respectivo cluster e líder
								 // Alterar para minimização
		
		double[] IToUpdate = new double[Run.DIMENSIONS]; 
		
		for(int i=0;i<IToUpdate.length;i++){
			double sum1 = 0.0;
			double sum2 = 0.0;
			IToUpdate[i] = 0.0;
			for (MOFish fish: this.fishCluster){
				if ( Double.compare(fish.lastAggregatedWeight, fish.aggregatedWeight) > 0) {
					sum1 += (fish.position[i] - fish.lastPosition[i]) * (fish.lastAggregatedWeight - fish.aggregatedWeight);
					sum2 += fish.lastAggregatedWeight - fish.aggregatedWeight;
				}	
			}
			if(sum2 != 0){
				IToUpdate[i] = sum1/sum2;
			}
			
		}
			
		I = IToUpdate;
		
	}

	public void moveFishColInstinctive() {
		if (!this.leader){
			moveFish(I);
		}
		
		
	}
	
	public void calculateB(){   // Recebe um peixe, seu respectivo cluster e líder

		double[] BToUpdate = new double[Run.DIMENSIONS]; 
		for(int i=0;i<BToUpdate.length;i++){
			BToUpdate[i] = 0.0;
			double sum1 = 0;
			double sum2 = 0;
			double sum3 = 0;
			for (MOFish fish: this.fishCluster){
				
				sum1 += fish.position[i] * (1/(fish.aggregatedWeight+1));
				sum2 += (1/(fish.aggregatedWeight+1));
				sum3 += (1/(fish.lastAggregatedWeight+1));
					
			}
			
			if(sum2 != 0){
				BToUpdate[i] = sum1/sum2;
			}
			
			if (Double.compare(sum2, sum3) > 0){
				sense = -1;
			} else {
				sense = 1;
			}
		}
			
		B = BToUpdate;
	}

	public void moveFishVolitive(double stepVol) {		//Já alterado para minimização
				
		double[] displacement = new double[Run.DIMENSIONS];
		double randNum;
		double distance = euclidianDistance(B, position);
		double constant;
		
		if (!this.leader) {
			if (sense < 0) {
				for(int i=0; i<displacement.length; i++){
					randNum = Run.generator.nextDouble();
					if (distance < 0.000001) {
						constant = 0;
					} else{
						constant = -stepVol * randNum/distance;
					}
					displacement[i] = (position[i] - B[i]) * constant;
				}
				moveFish(displacement);			
			} else {
				
				for(int i=0; i<displacement.length; i++){
					randNum = Run.generator.nextDouble();
					if (distance < 0.000001) {
						constant = 0;
					} else{
						constant = stepVol * randNum/distance;
					}
					displacement[i] = (position[i] - B[i]) * constant;
				}
				moveFish(displacement);
			}
		}
	}
	
	public double euclidianDistance(double[] a, double[] b){
		
		double total = 0;
		
		for(int i=0;i<a.length;i++){
			total += Math.pow((a[i] - b[i]), 2);		
		}
		
		total = Math.sqrt(total);
				
		return Math.max(total, 0.000001);
	}
		
	public void setRefPointIndex(int refPointIndexPassed) {
		refPointIndex = refPointIndexPassed;
		
	}

	public int getRefPointIndex() {
		// TODO Auto-generated method stub
		return this.refPointIndex;
	}
	
	public static void printList(ArrayList<double[]> listToPrint){
		System.out.println("Lista = " + listToPrint);
		for(double[] point: listToPrint){
			for(int i=0; i< point.length; i++) {
			System.out.println("Component[" + i + "] = " + point[i]);
			}
		}
	}
	
	public static void printListOfFishes(ArrayList<MOFish> listToPrint){
		System.out.println("Lista = " + listToPrint);
		for(MOFish fish: listToPrint){
			for(int i=0; i< fish.position.length; i++) {
			System.out.println("Component[" + i + "] = " + fish.position[i]);
			}
		}
	}
	
	public static void printArray(double[] arrayToPrint){
		System.out.println("Array = " + arrayToPrint);
			for(int i=0; i< arrayToPrint.length; i++) {
			System.out.println("Component[" + i + "] = " + arrayToPrint[i]);
			}
	}

	public void feedFishNotNormalizing() {
		
		this.updateWeight(this.fitness.clone());
		this.lastAggregatedWeight = this.aggregatedWeight;
		this.aggregatedWeight = this.aggregate(weight);
		
	}

	public void resetLeader() {

		this.leader = false;
		
	}

	public void setLeader() {

		this.leader = true;
		
	}
	
	public boolean getIsLeader() {

		return this.leader;
		
	}
}
