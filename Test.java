package wmoFSS;

import General.Run;

public class Test {
    public static void main(String [] args) {
    	
    	int GRID1 = 12;
    	int GRID2 = 1;
    	int NUMBER_OF_OBJECTIVES = 3;
    	
		int n1 = GRID1 + NUMBER_OF_OBJECTIVES - 1;
		int r = NUMBER_OF_OBJECTIVES - 1;
		double N1 = getFactorial(n1)/(getFactorial(n1 - r) * getFactorial(r));
		int n2 = GRID2 + NUMBER_OF_OBJECTIVES - 1;
		double N2 = getFactorial(n2)/(getFactorial(n2 - r) * getFactorial(r));
		double NUMBER_OF_REFPOINTS = N1 + N2;
		System.out.println("Fatorial n1 " + getFactorial(n1));
		System.out.println("Fatorial r " + getFactorial(r));
		System.out.println("Fatorial n1-r " + getFactorial(n1 - r));
		System.out.println("N1 " + N1);
		
		
		System.out.println("Fatorial n2 " + getFactorial(n2));
		System.out.println("Fatorial r " + getFactorial(r));
		System.out.println("Fatorial n2-r " + getFactorial(n2 - r));
		System.out.println("N2 " + N2);
		System.out.println("Number of Ref Points " + NUMBER_OF_REFPOINTS);
		

    }

    private static double getFactorial(int n) {
    	if (n < 0) {
    		System.out.println("NUMERO NEGATIVO!");
    	}
    	double factorial = 1.0;
    	for(double i=n; i>=1; i--){
    		factorial *= i;
    	}
	
    	return factorial;
	
    }
}
    	
    	
    	


