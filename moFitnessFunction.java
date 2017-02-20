package wmoFSS;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

import General.Run;

public class moFitnessFunction {
	public int NUMBER_OF_REFPOINTS;
	public ArrayList<double[]> referencePointsOutter;
	public ArrayList<double[]> referencePointsInner;
	public ArrayList<double[]> referencePoints;
	public ArrayList<double[]> paretoOptimalFront;
	//public double[] referencePointHV;
	
	public moFitnessFunction() throws FileNotFoundException {
		
		double n1 = Run.GRID1 + Run.NUMBER_OF_OBJECTIVES - 1;
		double r = Run.NUMBER_OF_OBJECTIVES - 1;
		double N1 = getFactorial(n1)/(getFactorial(n1 - r) * getFactorial(r));
		double n2 = Run.GRID2 + Run.NUMBER_OF_OBJECTIVES - 1;
		double N2 = 0;
		
		if (Run.GRID2 != 0) {
			N2 = getFactorial(n2)/(getFactorial(n2 - r) * getFactorial(r));
		}
		
		
		NUMBER_OF_REFPOINTS = (int) (N1 + N2);
		System.out.println("Number of Ref Points " + NUMBER_OF_REFPOINTS);
		
		referencePoints = new ArrayList<double[]>();
		
		referencePointsOutter = generatePoints(Run.GRID1);
		//printList(referencePointsOutter);
		
		referencePointsInner = generatePoints(Run.GRID2);
		//printList(referencePointsInner);
		
		referencePoints.addAll(referencePointsOutter);
		referencePoints.addAll(referencePointsInner);
		
		paretoOptimalFront = new ArrayList<double[]>();
		paretoOptimalFront = readParetoOptimalFront();
		//printList(paretoOptimalFront);
		
		//referencePointHV = new double[Run.NUMBER_OF_OBJECTIVES];
		
	
	}
	
	
	private double getFactorial(double n) {
		if (n < 0) {
			System.out.println("NUMERO NEGATIVO!");
		}
		double factorial = 1.0;
		for(double i=n; i>=1; i--){
			factorial *= i;
        }
		
		return factorial;
	}
	
	public ArrayList<double[]> getReferencePoints(){
		return referencePoints;
	}
	
	public double[] calculateFitness(double[] position) {
		int numberOfObjectives = Run.NUMBER_OF_OBJECTIVES;
		int numberOfVariables = Run.DIMENSIONS;
		
		double[] mofitness = new double[Run.NUMBER_OF_OBJECTIVES];

		double[] f = new double[numberOfObjectives];
		double[] x = new double[numberOfVariables] ;

		int k = numberOfVariables - numberOfObjectives + 1;

		for (int i = 0; i < numberOfVariables; i++) {
			x[i] = position[i] ;
		}
		
		switch(Run.TESTSET) {
		
		case 1: // DTLZ1
		// Cálculo do funcional g
		    double g = 0.0;
		    for (int i = numberOfVariables - k; i < numberOfVariables; i++) {
		      g += (x[i] - 0.5) * (x[i] - 0.5) - Math.cos(20.0 * Math.PI * (x[i] - 0.5));
		    }

		    g = 100 * (k + g);
		    
		    // Cálculo dos objetivos f1(x), f2(x), ...
		    for (int i = 0; i < numberOfObjectives; i++) {
		      f[i] = (1.0 + g) * 0.5;
		    }

		    for (int i = 0; i < numberOfObjectives; i++) {
		      for (int j = 0; j < numberOfObjectives - (i + 1); j++) {
		        f[i] *= x[j];
		      }
		      if (i != 0) {
		        int aux = numberOfObjectives - (i + 1);
		        f[i] *= 1 - x[aux];
		      }
		    }

		    for (int i = 0; i < numberOfObjectives; i++) {
		      mofitness[i] = f[i];
		    }
		    break;
		  
		case 2: //DTLZ 2   
		
			g = 0.0;
		    for (int i = numberOfVariables - k; i < numberOfVariables; i++) {
		      g += (x[i] - 0.5) * (x[i] - 0.5);
		    }

		    for (int i = 0; i < numberOfObjectives; i++) {
		      f[i] = 1.0 + g;
		    }

		    for (int i = 0; i < numberOfObjectives; i++) {
		      for (int j = 0; j < numberOfObjectives - (i + 1); j++) {
		        f[i] *= Math.cos(x[j] * 0.5 * Math.PI);
		      }
		      if (i != 0) {
		        int aux = numberOfObjectives - (i + 1);
		        f[i] *= Math.sin(x[aux] * 0.5 * Math.PI);
		      }
		    }

		    for (int i = 0; i < numberOfObjectives; i++) {
		    	mofitness[i] = f[i];
		    }
		    
		    break;
		
		case 3: //DTLZ 3
			g = 0.0;
		    for (int i = numberOfVariables - k; i < numberOfVariables; i++) {
		      g += (x[i] - 0.5) * (x[i] - 0.5) - Math.cos(20.0 * Math.PI * (x[i] - 0.5));
		    }

		    g = 100.0 * (k + g);
		    for (int i = 0; i < numberOfObjectives; i++) {
		      f[i] = 1.0 + g;
		    }

		    for (int i = 0; i < numberOfObjectives; i++) {
		      for (int j = 0; j < numberOfObjectives - (i + 1); j++) {
		        f[i] *= java.lang.Math.cos(x[j] * 0.5 * java.lang.Math.PI);
		      }
		      if (i != 0) {
		        int aux = numberOfObjectives - (i + 1);
		        f[i] *= java.lang.Math.sin(x[aux] * 0.5 * java.lang.Math.PI);
		      }
		    }

		    for (int i = 0; i < numberOfObjectives; i++) {
		      mofitness[i] = f[i];
		    }
		    
		    break;
		    
		case 4: //DTLZ 4
			double alpha = 100;
			g = 0.0;
		    for (int i = numberOfVariables - k; i < numberOfVariables; i++) {
		      g += (x[i] - 0.5) * (x[i] - 0.5);
		    }

		    for (int i = 0; i < numberOfObjectives; i++) {
		      f[i] = 1.0 + g;
		    }

		    for (int i = 0; i < numberOfObjectives; i++) {
		      for (int j = 0; j < numberOfObjectives - (i + 1); j++) {
		        f[i] *= java.lang.Math.cos(java.lang.Math.pow(x[j], alpha) * (java.lang.Math.PI / 2.0));
		      }
		      if (i != 0) {
		        int aux = numberOfObjectives - (i + 1);
		        f[i] *= java.lang.Math.sin(java.lang.Math.pow(x[aux], alpha) * (java.lang.Math.PI / 2.0));
		      }
		    }

		    for (int i = 0; i < numberOfObjectives; i++) {
		      mofitness[i] = f[i];
		    }
		    
		    break;
		    
		case 7: //DTLZ 7
			g = 0.0;
		    for (int i = numberOfVariables - k; i < numberOfVariables; i++) {
		      g += x[i];
		    }

		    g = 1 + (9.0 * g) / k;

		    System.arraycopy(x, 0, f, 0, numberOfObjectives - 1);

		    double h = 0.0;
		    for (int i = 0; i < numberOfObjectives - 1; i++) {
		      h += (f[i] / (1.0 + g)) * (1 + Math.sin(3.0 * Math.PI * f[i]));
		    }

		    h = numberOfObjectives - h;

		    f[numberOfObjectives - 1] = (1 + g) * h;

		    for (int i = 0; i < numberOfObjectives; i++) {
		      mofitness[i] = f[i];
		    }
		    
		    break;
		}
		return mofitness;
	}
	
	
	/**
	 * Generates the reference points (weights) for the given number of
	 * divisions.
	 * 
	 * @param divisions the number of divisions
	 * @return the list of reference points
	 */
	private ArrayList<double[]> generatePoints(int divisions) {
		ArrayList<double[]> result = new ArrayList<double[]>();
		double[] weight = new double[Run.NUMBER_OF_OBJECTIVES];
		
		if (divisions != 0) {
			generateRecursive(result, weight, Run.NUMBER_OF_OBJECTIVES, divisions, divisions, 0);
		}
		
		return result;
	}
	
	/**
	 * Generate reference points (weights) recursively.
	 * 
	 * @param weights list storing the generated reference points
	 * @param weight the partial reference point being recursively generated
	 * @param numberOfObjectives the number of objectives
	 * @param left the number of remaining divisions
	 * @param total the total number of divisions
	 * @param index the current index being generated
	 */
	private void generateRecursive(ArrayList<double[]> weights,
			double[] weight, int numberOfObjectives, int left, int total, int index) {
		if (index == (numberOfObjectives - 1)) {
			weight[index] = (double)left/total;
			weights.add(weight.clone());
		} else {
			for (int i = 0; i <= left; i += 1) {
				weight[index] = (double) i / total;
				generateRecursive(weights, weight, numberOfObjectives, left - i, total, index + 1);
			}
		}
	}
	
	private ArrayList<double[]> readParetoOptimalFront() throws FileNotFoundException {
    	String file = new String();
    	//file = "DTLZ";
    	//file = file.concat(Integer.toString(Run.TESTSET));
    	//file = file.concat(Integer.toString(Run.NUMBER_OF_OBJECTIVES));
    	file = Run.file.concat(".txt");
    	System.out.println(file);
    	Scanner sc = new Scanner(new File(file));
    	sc.useDelimiter("\n");
    	ArrayList<double[]> paretoOptimalFront = new ArrayList<double[]>();

    	while (sc.hasNextLine()) {
    	    String[] parts = sc.nextLine().split(" "); // split each line by " "
    	    double[] temp = new double[Run.NUMBER_OF_OBJECTIVES];
    	    int i = 0;
    	    for (String s : parts) {
    	    	temp[i] = (double) Double.parseDouble(s);
    	    	i += 1;
    	    }
    	    paretoOptimalFront.add(temp);
    	}
    	sc.close();
    	return paretoOptimalFront;
	}
	
	
	/*public double[] getReferencePointHV(){
		switch (Run.TESTSET) {
		case 1:
			for (int i = 0; i < Run.NUMBER_OF_OBJECTIVES; i++) {
				referencePointHV[i] = 1.0;
			}
			break;
		case 2:
			for (int i = 0; i < Run.NUMBER_OF_OBJECTIVES; i++) {
				referencePointHV[i] = 2.0;
			}
			break;
		case 3:
			for (int i = 0; i < Run.NUMBER_OF_OBJECTIVES; i++) {
				referencePointHV[i] = 2.0;
			}
			break;
		case 4:
			for (int i = 0; i < Run.NUMBER_OF_OBJECTIVES; i++) {
				referencePointHV[i] = 2.0;
			}
			break;
		case 7:
			for (int i = 0; i < Run.NUMBER_OF_OBJECTIVES; i++) {
				referencePointHV[i] = 2.0;
			}
			break;	
		}	
		return referencePointHV;
	}*/
	
	public static void printList(ArrayList<double[]> listToPrint){
		System.out.println("Lista = " + listToPrint);
		for(double[] point: listToPrint){
			for(int i=0; i< point.length; i++) {
			System.out.println("Component[" + i + "] = " + point[i]);
			}
		}
	}
}


