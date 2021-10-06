

package jmetal.problems;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

import java.util.*;
import java.lang.Double.*;
import java.io.*;
import java.util.ArrayList;  
import java.util.Collections; 


//import flanagan.integration.Integration;
//import flanagan.integration.IntegralFunction;
 
public class ReliabilityProblem97BaseIN extends Problem  {    
    
 private static final double Infinty = 0;
 public static int indexIterationFC=0; 
/**
  * Constructor.
  * Creates a default instance of the ReliabilityProblem problem
  * @param solutionType The solution type must "Real" or "BinaryReal".
  */
	
//number of debuggers	
int m = 8;

//number of functionality 
int K = 8;

//istanza nominale
double B = 2500.0;
//istanza 0 
//double B = 2927;
//double B = 2796;
//double 2736;
//double B = 2815;
//tolerance level for the estimation error e of Monte-Carlo 
double Tol = 0.5;

double A = 0.8; 
double alpha= 0.5; 
double kappa= 0.05;


//valore della media e varianza delle distribuzioni dai paremetri dell'SRGM delle funzionalita' 
//per ogni parametro dovremmo definire media e varianza
//assumiamo P parametri
//int P = 2;
double [] weight ={0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
//double [] weight ={0.7, 0.2, 0.1, 0, 0, 0, 0, 0};
double Rmin = 0.97;
double F = 1 - Rmin;

//queste sono le variabili booleane del modello. Sono pari ad uno se il deb e' coinvolto per la funz. k 
int [][] x = new int [m][K];

//ESEMPIO di EXPONENZIALE
//range per uncertainty

//double [][] range = {{0.2475, 0.2676}, {0.211, 0.2166}, {0.1535, 0.1656}, {0.2871, 0.3288}, {0.5665, 0.7159}, {0.1971, 0.22511}, {0.1828, 0.2231}, {0.2178, 0.2661}}; 
// ISTANZA NOMINALE
double [][] range = {{0.02475, 0.02676}, {0.0211, 0.02166}, {0.01535, 0.01656}, {0.02871, 0.03288}, {0.05665, 0.07159}, {0.01971, 0.022511}, {0.01828, 0.02231}, {0.02178, 0.02661}}; 

//allargato di 0.01
//double [][] range = {{0.01475, 0.3676}, {0.0111, 0.3166}, {0.00535, 0.2656}, {0.01871, 0.4288}, {0.04665, 0.8159}, {0.00971, 0.12511}, {0.00828, 0.3231}, {0.01178, 0.3661}};
//media numero di ore per lavorazione dei bug

//istanza nominale
double [] delta = {0, 0, 0,0, 0,0, 0, 0};
//ISTANCE 0
//double [] delta = {2.134829828250902, 1.8140908580842985, 1.8209372027914874, 1.9470523571306995, 1.0149727390964, 1.088805084972804, 2.1246697319322005, 1.006704130364219};
//ISTANCE 1
//double [] delta = {1.9920736340041556, 2.061382395392698, 1.9416365045000017, 2.0968452415984493, 1.0656857494860985, 1.0955427418019002, 2.091474979098764, 1.09192986712266}; 
//Istance 2
//double [] delta = {2.009890604076438, 2.084913544283662, 1.809841980423505, 2.1496167738008185, 0.9636328339332502, 0.9082543197642979, 2.0911179701541394, 1.0379026797685593}; 
//Instance 3
//double [] delta = {1.8108457962996947, 1.94296447344517, 1.86072665136233, 1.9954944960142837, 0.9122909717732326, 0.9938360225857581, 1.8631986622455123, 0.9848110266713678};
//numero medio di ore al giorno che un debugger puo' lavorare sulla funzionalita' k 
//istanza nominale

//ISTANZA 0

double [][] mTh = {{0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24},
		{0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24},
		{0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24},
		{0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24},
		{0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24},
		{0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24},
		{0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24},
		{0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24},
		{0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24, 0/24},
};


/*double [][] mTh = {{4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24},
		{4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24},
		{4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24},
		{4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24},
		{4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24},
		{4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24},
		{4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24},
		{4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24},
		{4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24, 4/24}
		
};
*/
//numero atteso iniziale di difetti
//Istanza Nominale
double [] omega = {204.59, 215.39, 60.60, 39.17, 33.74, 21.31, 20.42, 22.41}; 
//Istanza 0
//double [] omega = {187.49422598372635, 224.96981657550785, 62.51480251069716, 36.51764556119514, 30.89970338547098, 20.144836793229874, 20.745915530588654, 21.21114579589483};  
//Istanza 1
//double [] omega = {184.86456163820338, 216.2959789862524, 65.0898316846495, 36.26481091544923, 33.40123178084089, 20.994923756708726, 20.35511962639653, 23.167463961412068}; 
//Istanza 2
//double [] omega = {192.85591586453987, 229.66988744880894, 57.1505113721635, 35.794652073760886, 36.02851287346513, 19.392113220718215, 18.869498250968302, 24.170874372449305}; 
//Istanza 3
//double [] omega = {220.53338605825903, 209.21012855212464, 65.59988562508767, 39.20758535311831, 36.07520672808535, 21.25408341442324, 18.934180984142138, 22.47815616296878}; 


//cost of a personday
double costmd = 1;

public void CreazioneIstanze(PrintStream p) {

	int i, k, j;
	double app, app1, app2, Nist=30;
	
	
	//calcolo il 10 per cent di delta
	
for (k=0; k< Nist; k++) {
	
	p.println("\n");
	p.println("ISTANZA" + k);
/*	
	p.println("valori di delta");
	for (i=0; i<m; i++) {
	
		app = (delta[i]*10)/100;
		
		//genero un numero tra 1 e 10. Se il num e' maggiore di 5 incremento, altrimenti decremento
		
		app1 = ( (int)(Math.random() * 10));
		
		if (app1 > 5) app2 = delta[i] + (app * ( (double)(Math.random() % (double) 1)) );
		else app2 = delta[i] - (app * ( (double)(Math.random() % (double) 1)) );
		
	     p.print(+app2+", ");
		    
	}
	*/	
		p.println("\n");
		p.println("valori di omega");
	
		
		for (i=0; i<K; i++) {
		
			app = (omega[i] * 10)/100;
			
			//genero un numero tra 1 e 10. Se il num e' maggiore di 5 incremento, altrimenti decremento
			
			app1 = ( (int)(Math.random() * 10));
			
			if (app1 > 5) app2 = omega[i] + (app * ( (double)(Math.random() % (double) 1)) );
			else app2 = omega[i] - (app * ( (double)(Math.random() % (double) 1)) );
			
			p.print(+app2+", ");
			    
		}//end for per omega
		
		
		p.println("\n");
		p.println("valori di B");
		
		app = (B * 20)/100;
		
		//genero un numero tra 1 e 10. Se il num e' maggiore di 5 incremento, altrimenti decremento
		
		app1 = ( (int)(Math.random() * 10));
		
		if (app1 > 5) app2 = B + (app * ( (double)(Math.random() % (double) 1)) );
		else app2 = B - (app * ( (double)(Math.random() % (double) 1)) );
		
		p.println(app2);
	
	
	}//end for sulle istanze
	
}



public class Randoms extends java.util.Random {

	  public Randoms (int seed) {
	    super(seed);
	  }
	  
	  public Randoms () {
	    super();
	  }
	  
	  
	  /** Return random integer from Poission with parameter lambda.  
	   * The mean of this distribution is lambda.  The variance is lambda. */
	  public synchronized int nextPoisson(double lambda) {
	    int i,j,v=-1;
	    double l=Math.exp(-lambda),p;
	    p=1.0;
	    while (p>=l) {
	      p*=nextUniform();
	      v++;
	    }
	    return v;
	  }

	  /** Return nextPoisson(1). */
	  public synchronized int nextPoisson() {
	    return nextPoisson(1);
	  }

	  /** Return a random boolean, equally likely to be true or false. */
	  public synchronized boolean nextBoolean() {
	    return (next(32) & 1 << 15) != 0;
	  }

	  /** Return a random boolean, with probability p of being true. */
	  public synchronized boolean nextBoolean(double p) {
	    double u=nextUniform();
	    if(u < p) return true;
	    return false;
	  }

	  /** Return a random BitSet with "size" bits, each having probability p of being true. */
	  public synchronized BitSet nextBitSet (int size, double p)
	  {
	    BitSet bs = new BitSet (size);
	    for (int i = 0; i < size; i++)
	      if (nextBoolean (p)) {
	        bs.set (i);
	      }
	    return bs;
	  }

	  /** Return a random double in the range 0 to 1, inclusive, uniformly sampled from that range. 
	   * The mean of this distribution is 0.5.  The variance is 1/12. */
	  public synchronized double nextUniform() {
		//qua ho usato nextInt invece di next()
		  return Math.random();
	//	  long l = ((long)(nextInt(26)) << 27) + nextInt(27);
	  //  return l / (double)(1L << 53);
	  }

	  /** Return a random double in the range a to b, inclusive, uniformly sampled from that range.
	   * The mean of this distribution is (b-a)/2.  The variance is (b-a)^2/12 */
	  public synchronized double nextUniform(double a,double b) {
	    return a + (b-a)*nextUniform();
	  }

	  /** Draw a single sample from multinomial "a". */
	  public synchronized int nextDiscrete (double[] a) {
	    double b = 0, r = nextUniform();
	    for (int i = 0; i < a.length; i++) {
	      b += a[i];
	      if (b > r) {
	        return i;
	      }
	    }
	    return a.length-1;
	  }

	  /** draw a single sample from (unnormalized) multinomial "a", with normalizing factor "sum". */
	  public synchronized int nextDiscrete (double[] a, double sum) {
	    double b = 0, r = nextUniform() * sum;
	    for (int i = 0; i < a.length; i++) {
	      b += a[i];
	      if (b > r) {
	        return i;
	      }
	    }
	    return a.length-1;
	  }

	  private double nextGaussian;
	  private boolean haveNextGaussian = false;

	  /** Return a random double drawn from a Gaussian distribution with mean 0 and variance 1. */
	  public synchronized double nextGaussian() {
	    if (!haveNextGaussian) {
	      double v1=nextUniform(),v2=nextUniform();
	      double x1,x2;
	      x1=Math.sqrt(-2*Math.log(v1))*Math.cos(2*Math.PI*v2);
	      x2=Math.sqrt(-2*Math.log(v1))*Math.sin(2*Math.PI*v2);
	      nextGaussian=x2;
	      haveNextGaussian=true;
	      return x1;
	    }
	    else {
	      haveNextGaussian=false;
	      return nextGaussian;
	    }
	  }

	  /** Return a random double drawn from a Gaussian distribution with mean m and variance s2. */
	  public synchronized double nextGaussian(double m,double s2) {
	    return nextGaussian()*Math.sqrt(s2)+m;
	  }

	  // generate Gamma(1,1)
	  // E(X)=1 ; Var(X)=1
	  /** Return a random double drawn from a Gamma distribution with mean 1.0 and variance 1.0. */
	  public synchronized double nextGamma() {
	    return nextGamma(1,1,0);
	  }

	  /** Return a random double drawn from a Gamma distribution with mean alpha and variance 1.0. */
	  public synchronized double nextGamma(double alpha) {
	    return nextGamma(alpha,1,0);
	  }

	  /* Return a sample from the Gamma distribution, with parameter IA */
	  /* From Numerical "Recipes in C", page 292 */
	  public synchronized double oldNextGamma (int ia)
	  {
	    int j;
	    double am, e, s, v1, v2, x, y;

	    assert (ia >= 1) ;
	    if (ia < 6) 
	    {
	      x = 1.0;
	      for (j = 1; j <= ia; j++)
	        x *= nextUniform ();
	      x = - Math.log (x);
	    }
	    else
	    {
	      do
	      {
	        do
	        {
	          do
	          {
	            v1 = 2.0 * nextUniform () - 1.0;
	            v2 = 2.0 * nextUniform () - 1.0;
	          }
	          while (v1 * v1 + v2 * v2 > 1.0);
	          y = v2 / v1;
	          am = ia - 1;
	          s = Math.sqrt (2.0 * am + 1.0);
	          x = s * y + am;
	        }
	        while (x <= 0.0);
	        e = (1.0 + y * y) * Math.exp (am * Math.log (x/am) - s * y);
	      }
	      while (nextUniform () > e);
	    }
	    return x;
	  }

	  
	  /** Return a random double drawn from a Gamma distribution with mean alpha*beta and variance alpha*beta^2. */
	  public synchronized double nextGamma(double alpha, double beta) {
	    return nextGamma(alpha,beta,0);
	  }

	  /** Return a random double drawn from a Gamma distribution with mean alpha*beta+lamba and variance alpha*beta^2. */
	  public synchronized double nextGamma(double alpha,double beta,double lambda) {
	    double gamma=0;
	    if (alpha <= 0 || beta <= 0) {
	      throw new IllegalArgumentException ("alpha and beta must be strictly positive.");
	    }
	    if (alpha < 1) {
	      double b,p;
	      boolean flag=false;
	      b=1+alpha*Math.exp(-1);
	      while(!flag) {
	        p=b*nextUniform();
	        if (p>1) {
	          gamma=-Math.log((b-p)/alpha);
	          if (nextUniform()<=Math.pow(gamma,alpha-1)) flag=true;
	        }
	        else {
	          gamma=Math.pow(p,1/alpha);
	          if (nextUniform()<=Math.exp(-gamma)) flag=true;
	        }
	      }
	    }
	    else if (alpha == 1) {
	      gamma = -Math.log (nextUniform ());
	    } else {
	      double y = -Math.log (nextUniform ());
	      while (nextUniform () > Math.pow (y * Math.exp (1 - y), alpha - 1))
	        y = -Math.log (nextUniform ());
	      gamma = alpha * y;
	    }
	    return beta*gamma+lambda;
	  }

	  /** Return a random double drawn from an Exponential distribution with mean 1 and variance 1. */
	  public synchronized double nextExp() {
	    return nextGamma(1,1,0);
	  }

	  /** Return a random double drawn from an Exponential distribution with mean beta and variance beta^2. */
	  public synchronized double nextExp(double beta) {
	    return nextGamma(1,beta,0);
	  }

	  /** Return a random double drawn from an Exponential distribution with mean beta+lambda and variance beta^2. */
	  public synchronized double nextExp(double beta,double lambda) {
	    return nextGamma(1,beta,lambda);
	  }

	  /** Return a random double drawn from an Chi-squarted distribution with mean 1 and variance 2. 
	   * Equivalent to nextChiSq(1) */
	  public synchronized double nextChiSq() {
	    return nextGamma(0.5,2,0);
	  }

	  /** Return a random double drawn from an Chi-squared distribution with mean df and variance 2*df.  */
	  public synchronized double nextChiSq(int df) {
	    return nextGamma(0.5*(double)df,2,0);
	  }

	  /** Return a random double drawn from an Chi-squared distribution with mean df+lambda and variance 2*df.  */
	  public synchronized double nextChiSq(int df,double lambda) {
	    return nextGamma(0.5*(double)df,2,lambda);
	  }

	  /** Return a random double drawn from a Beta distribution with mean a/(a+b) and variance ab/((a+b+1)(a+b)^2).  */
	  public synchronized double nextBeta(double alpha,double beta) {
	    if (alpha <= 0 || beta <= 0) {
	      throw new IllegalArgumentException ("alpha and beta must be strictly positive.");
	    }
	    if (alpha == 1 && beta == 1) {
	      return nextUniform ();
	    } else if (alpha >= 1 && beta >= 1) {
	      double A = alpha - 1,
	              B = beta - 1,
	              C = A + B,
	              L = C * Math.log (C),
	              mu = A / C,
	              sigma = 0.5 / Math.sqrt (C);
	      double y = nextGaussian (), x = sigma * y + mu;
	      while (x < 0 || x > 1) {
	        y = nextGaussian ();
	        x = sigma * y + mu;
	      }
	      double u = nextUniform ();
	      while (Math.log (u) >= A * Math.log (x / A) + B * Math.log ((1 - x) / B) + L + 0.5 * y * y) {
	        y = nextGaussian ();
	        x = sigma * y + mu;
	        while (x < 0 || x > 1) {
	          y = nextGaussian ();
	          x = sigma * y + mu;
	        }
	        u = nextUniform ();
	      }
	      return x;
	    } else {
	      double v1 = Math.pow (nextUniform (), 1 / alpha),
	              v2 = Math.pow (nextUniform (), 1 / beta);
	      while (v1 + v2 > 1) {
	        v1 = Math.pow (nextUniform (), 1 / alpha);
	        v2 = Math.pow (nextUniform (), 1 / beta);
	      }
	      return v1 / (v1 + v2);
	    }
	  }
	  
}  

//SECONDA VERSIONE PER IL CAMPIONAMENTO DA DISTRIBUZIONE ESPONENZIALE

public class ExponentialRV {
	
	private double mu= 1.;

    public ExponentialRV (double mu){
	             	this.mu = mu;     
    	}

    public double nextDouble(){
    	 Random r = new Random ();
    	return  - mu * Math.log(r.nextDouble());     
 }

}


public double MeanConfidenceInterval(boolean low, double me, double var) {

	 double stddev = Math.sqrt(var);
	 
	 //lower bound
	 if (low == true)		 
		//1.96 e' legato al 95%
		  return  me - 1.96 * stddev;
	
    //upper bound
	 else return me + 1.96 * stddev;

	}

public class fuctionIntegral {
	
	    private double mu1, beta1, Yk1, Tk1; 
	    int k1;

	    public void setMU(double mu){
    	  	this.mu1 = mu;     
      }
 
      public void setBeta(double beta){
            this.beta1 = beta;        
      }
      
      public void setYT(double Y){
          this.Yk1 = Y;        
    }
      
      
      public void setT(double T){
          this.Tk1 = T;        
    }
       
      public void setK(int k){
          this.k1 = k;        
    }  
      
	    public double power(double number,int elevation){
			double ret=1;
			for(int i=elevation;i>0;i--){
				ret *= number;
				}
				return ret;
			}
        public double f(double x)                                       
        {      double y=0;
              // y = mu1 * Math.exp(mu1*x) * getMeanValueFunction(x,shape1,scale1); 
        //questa e' la funzione che sta dentro l'integrale di m_ck
        
         y= omega[k1]* mu1 * Math.exp(mu1*x) * (1 - Math.exp(-(beta1 * Yk1)));
         
       //  System.out.println("y vale "+ y);
        //System.out.println("x vale "+ x);
       // y=power(x,2);
                return y;
        }
         public double trapArea(double b1, double b2, double h)         
        {
                return (b1+b2)*h/2;
        }
 
        public double trapRule(double x0, double x1, int div)          
        {
                double area = 0;
                double h = (x1 - x0) / div;                                                              
                for (int i=0; i<div-1; i++)                            
                {                                                    
                        area += trapArea(f(x0), f(x0+h), h);            
                        x0 += h;
                }
                area += trapArea(f(x0), f(x1), x1-x0);
                return area;
 
        }
}

//nel caso esponenziale

public double mdk(int k, double betak, double Yk) {

	return omega[k]*(1 - Math.exp(-(betak * Yk)));
}


public double mck(int k, double betak, double muk, double Yk, double Tk) {

fuctionIntegral  integral = new fuctionIntegral();

integral.setMU(muk);
integral.setBeta(betak);
integral.setYT(Yk);
integral.setT(Tk);
integral.setK(k);

double number = integral.trapRule(0,Tk,5000);


//System.out.println("number vale "+ number);

return Math.exp(-(muk*Tk))* number;

}

/** Monte-Carlo simulation with dynamic stopping criteria */
// e is the error
//questa funz serve per calcolare il lower bound della prima funzione obiettivo
//calcola il lower bound della somma dei bug corretti per ogni funzionalita'

//SENZA I DEBUGGER
public double MonteCarloF1(boolean lconf,  double [][] C, double [] Y) {
	
	
	int i=0;
	int j=0;
	//finestra di campioni che consideriamo per andare a controllare l'errore
	int h = 10;
	int k;
	int per;
	double valper;
	int t;
	int d;
	
	double Tk = 0;
	 
	double aver=0.0D;
	double sumSquares=0.0D, stdev;
	
	double averTot=0, vrnTot= 0;
	int nRun=0;
	
	double e=0.6;
	
	double [] mu = new double [K];
	//nel caso esponenziale abbiamo solo un parametro, il beta
	double [] pSRGM = new double [K]; 

	//Z dobbiamo inizializzarla alla prima parte della formula dell'errore
	//1.96 e' legato al 95%
	double Z = 1.96;
	double totBug = 0.0;
	int flag =0;
	//caso di exponential
	double [] beta= new double [K];
	
	double app6, app7, app1;
	boolean check = true;
	
	//lista dei valori di reliability per ogni MC run
	//LinkedList<Double> L = new LinkedList<Double>();
	List<Double> L = new ArrayList <Double>();
	
	//lista del 20th percentile ottenuto per una finestra di h MC run 
	List<Double> Perc = new ArrayList <Double>();
    
	
	//campioniamo h valori alla volta, al h-simo andiamo a controllare l'errore
	
	while (e > 0.00001){
	//	while (nRun < 400){
		if (i < h){
				
			check = true;
			
			for (k = 0; k < K; k++){	
			//campioniamo i valori di Thkd
				//System.out.println("campionamento inizio PER FUNZIONALITA " + k); 
				mu[k]=0;
				//con la prima versione di campionamento
				Randoms r = new Randoms();
				
				
				//passiamo alla funzione che campiona nextUniform il range di ogni parametro
			    pSRGM [k] = r.nextUniform(range[k][0], range[k][1]);
			  //  if(k==0) System.out.println("range["+k+"],0: "+range[k][0]);
			   // if(k==0) System.out.println("range["+k+"],1: "+range[k][1]);
			    //if(k==0) System.out.println("psrgm["+k+"]: "+pSRGM [k]);
			    //caso di Exponential
			    beta[k]= pSRGM[k];
			    
			    // System.out.println(")prova BETA "+ beta[k]);
				//calcolo il mdk
			    app1 = mdk(k, beta[k],Y[k]);
			     
			     // System.out.println("app1 vale "+ app1); 
			    //campioniamo un valore del numero medio di ore di lavorazione di un bug con distrib esponenziale				
		//	    ExponentialRV Ex = new ExponentialRV(delta[k]);
			    app6=0;
		    
			   mu[k] = 0;
		 
				
			if (Y[k] > 0.0 && Y[k] < (B- 0.1)) 
	        	{
				
					Tk = (-1/(alpha*kappa))* Math.log((Math.pow((B/Y[k]),kappa) - 1.0)/A);        	      	
	        	
	        	}
			
	       else {if (Y[k] >= (B - 0.001)) Tk = (-1/(alpha*kappa))* Math.log((Math.pow((B/(B - 0.001)),kappa) - 1.0)/A);  
	       }
				
			//System.out.println("Tk vale "+ Tk);
				//CALCOLARE TK e PASSARLO IN MCK
			   app7= app1;//mck(k, beta[k], mu[k], Y[k], Tk);
			   
			   //System.out.println("mck vale "+ app7);
			    //calcolo il total bug corretti del sistema sommando quelli di ogni funzionalita'
			    //questo rappresenta il valore ottenuto con una run di MC
				if ((Double.isNaN(app7)) || (app7 == Infinty)) check = false; 
				else totBug = app7 + totBug;
		  	}
			
			if (check == true) {
			//mi porto dietro il valore della media e della varianza di tutte le run di MC
				
		//		System.out.println("\nTOT BUG \n"+totBug);
			averTot = averTot + totBug; 
        	vrnTot = vrnTot + (totBug * totBug);  // SOMMA DEI QUADRATI; sotto diventa varianza
						
        	//con i valori campionati calcoliamo il numero di bug rimanenti
		  	//memorizzo dentro la lista il valore
			//L.addLast(totBug);
		  	
        	L.add(totBug);
        	//ordino la lista
		    //System.out.println("lista ordinata ");		  	
		  	Collections.sort(L);
		  	//System.out.println("List value after sort: "+L); 
		  	
		    //calcoliamo il 20th percentile della lista
		  	//indice del 20th percentile    
//		 	per = (L.size()*20)/100;
		 	per = (int)Math.round(L.size()*0.05); // quinto percentile 
		 	valper= L.get(per);
		 	//memorizzo dentro la lista Perc il 20th percentile trovato
			Perc.add(valper);
			//System.out.println("percentile Lista ");
			//System.out.println("List PERC value: "+Perc);		  	
			totBug=0;
		   	i=i+1;
		   	nRun = nRun + 1;
		   	j=j+1;
		  }//end if check
		
		}//fine if i
	if (j >= h){
					   
		aver=0.0D;
		sumSquares=0.0D;
		for (t = Perc.size()- h; t < Perc.size(); t++){
	        		        		  			         
	        		 double aDouble = Perc.get(t);	        		 
	        	//	  System.out.println("aDouble Vale \n"+ aDouble);
	        		 aver = aver + aDouble; 
	        		 sumSquares = sumSquares + (aDouble * aDouble); 
	        		// System.out.println(")prova AVER "+ aver);
			    }
			if (aver != 0) {
		  		//calcoliamo l'errore
		  		//stdev = Math.sqrt((sumSquares- aver*aver/h) / (h - 1));
		  		stdev = Math.sqrt((sumSquares/h - (aver/h)*(aver/h)))* (double)h/(h - 1);
		  	//  System.out.println("aSTDEV  \n"+ stdev);
		  		if ((aver/h) != 0) e = 3.92 * (stdev/ ((Math.sqrt(h)* (aver/h))) );
		  	
		  	
		  
			}
		  		i=0;
	}//fine if j
		  
}//fine while
	
	/* STAMPA PER ERRORE
	System.out.println("Exited after nRUn: " +nRun);
	*/
	//restituisco l'intervallo di confidenza di tutte le MC run
	
	//System.out.println( "il lower bound inter conf"+ MeanConfidenceInterval(lconf,averTot/nRun,vrnTot/nRun));
	averTot = averTot/nRun;
	vrnTot = (vrnTot/nRun - Math.pow(averTot, 2))*((double)nRun/(nRun-1));
	
	//System.out.println("average "+averTot);
	//System.out.println("variance  "+vrnTot);
	//System.out.println("\n CONFIDENCE BUG "+MeanConfidenceInterval(lconf,averTot,vrnTot/nRun));
	if (nRun > 0) return MeanConfidenceInterval(lconf,averTot,vrnTot/nRun); // PASSO IL QUADRATO DELLO STANDARD ERROR
	
	else return 0;
    
}

//Monte Carlo simulazione per la funzione obiettivo 3 sul costo

public double MonteCarloF2(boolean lconf,  double [][] C, double [] Y) {
	
	int i=0;
	int j=0;
	//finestra di campioni che consideriamo per andare a controllare l'errore
	int h = 10;
	int k;
	int p;
	int t;
	int d;
	double aver=0.0D;
	double sumSquares=0.0D, stdev;
	double e=0.6D;
	double cost1=0.0D, app1=0.0D, app2=0.0D, app7;
	double averTot=0, vrnTot= 0;
	int per;
	double valper;
	int nRun = 0, checkval = 0;
	boolean check;
	double [] mu = new double [K];
	double [] pSRGM = new double [K]; 

	//caso di exponential
	double [] beta= new double [K];
	
	//Z dobbiamo inizializzarla alla prima parte della formula dell'errore
	//1.96 e' legato al 95%
	double Z = 1.96;
	double Tk = 0;
	
	//ATTENZIONE da stimare
	double C1 = 0, C2 = 80, C3 = 60; 
	double costTot = 0.0D, app=0.0, app6;

	//lista dei valori dei costi per ogni MC run
	//LinkedList<Double> L = new LinkedList<Double>();
	List<Double> L = new ArrayList <Double>();
		
	//lista del 20th percentile ottenuto per una finestra di h MC run 
	List<Double> Perc = new ArrayList <Double>();
	    		
	//campioniamo h valori alla volta, al h-simo andiamo a controllare l'errore
	
	while(e > 0.00001){
    
		if (i < h){
			
			check = true;
			
			for (k = 0; k < K; k++){
				
				app=0.0;
               //cost1=0.0;		
                
                mu[k]= 0;
				Randoms r = new Randoms();
				//campioniamo i valori di beta_k (cioe' del tasso di detection dell'SRGM)
				//caso di EXP usiamo uniforme
			    //passiamo alla funzione che campiona nextUniform il range di ogni parametro
			    pSRGM [k] = r.nextUniform(range[k][0], range[k][1]);			      		
			    //caso di Exponential
			    beta[k]= pSRGM[k];
			    
				//calcolo il mdk
			    app1 = mdk(k, beta[k],Y[k]);
			    
			    /*for (d = 0; d < m; d++){
					
			    	//calcolo la prima parte della funzione costo C1*_1k * mck(t) con eq 18			    	
					cost1 = (C[d][k] * app1 * costmd) + cost1;													
				}*/
				
			   // System.out.println("cost1 vale " + cost1);
			   //campioniamo un valore del numero medio di ore di lavorazione di un bug con distrib esponenziale				
				  ExponentialRV Ex = new ExponentialRV(delta[k]);
				  app6=Ex.nextDouble();
				  
					//System.out.println("app6 vale "+ app6);
					//app6=10;
				  
				  //vedere in caso chiudere
			
				  if (app6 < 0.5){
				   app6 = 1;
					}
			    
				//calcolo mu[k] con equazione 3 
				
				//for (d = 0; d < m; d++){	
					
						//mu[k]= (mTh[d][k]*Y[k]*C[d][k])/(app6*app1) + mu[k];
						//mu[k] = ( (((app6 * app1 * C[d][k])/ Y[k]) * mTh[d][k]) /app1) + mu[k];
						//mu[k] = ( (((app1 * C[d][k])/ Y[k]) * mTh[d][k]) /(app6* app1*C[d][k])) + mu[k];
						//if(C[d][k] > 0) mu[k] = ((mTh[d][k]/app6)/Y[k]) + mu[k];
						mu[k] = 0;//(1/app6);
						
			//}
				
				if (Y[k] > 0.0 && Y[k] < (B- 0.1)) 
			       
				{
						
				   Tk = (-1/(alpha*kappa))* Math.log((Math.pow((B/Y[k]),kappa) - 1.0)/A);        	      	
			        	
			     	}
					
			       else {if (Y[k] >= (B - 0.001)) Tk = (-1/(alpha*kappa))* Math.log((Math.pow((B/(B - 0.001)),kappa) - 1.0)/A);  
			       }
						
					//System.out.println("\nTk in f2 vale "+ Tk);
					//System.out.println(" Yk in f2 vale "+ Y[k]);
					//System.out.println(" delta in f2 vale "+ app6);
				    //CALCOLARE TK e PASSARLO IN MCK
				   app7= app1 ;//mck(k, beta[k], mu[k], Y[k], Tk);
				   
				   //System.out.println("mck in f2 vale "+ app7);
				   // System.out.println("omega in f2 vale "+ omega[k]);
					   
			       // C2= cost1 + 0.5;
				   if ((Double.isNaN(app7)) || (app7 == Infinty)) check = false; 
			
				   		else costTot= (C2 * (app6/24) * (omega[k] - app7) + C3 * Y[k]/24) + costTot;  	
			
				   	//System.out.println("il costo in f2 vale "+ costTot);	
		       }
			
			  if (check == true) {
				
				//mi porto dietro il valore della media e della varianza di tutte le run di MC
				averTot = averTot + costTot; 
				vrnTot = vrnTot + (costTot * costTot); 
			
				L.add(costTot);
				//ordino la lista
				//System.out.println("lista ordinata ");		  	
				Collections.sort(L);
				//System.out.println("List value affffffter sort: "+L); 
        	       	
				per = (int)Math.round(L.size()*0.05); // quinto percentile 
			 	valper= L.get(per);
			 	//memorizzo dentro la lista Perc il 20th percentile trovato
				Perc.add(valper);
				//System.out.println("percentile Lista ");
				//System.out.println("List PERC value: "+Perc);	
			
				costTot=0;
				i=i+1;
				
				nRun = nRun + 1;
			   	j=j+1;
			}//end if check 
		}//fine if i
	if (j >= h){
		
			//calcoliamo la media e varianza dei bug corretti
		//System.out.println("CALCOLO MEDIA E VARIANZA");
		  		aver=0;
		  		sumSquares=0;
		  		for (t = L.size()- h; t < L.size(); t++){
	        		        	
	        		        	double aDouble = Perc.get(t);
	        		        	aver = aver + aDouble; 
	        		        	sumSquares = sumSquares + (aDouble * aDouble); 
	        		        	//System.out.println(")prova AVER "+ aver);
			    			}
		  		//calcoliamo l'errore
		  		if (aver != 0) {
	 			//stdev = Math.sqrt((sumSquares- aver*aver/h) / (h - 1));
		  			stdev = Math.sqrt((sumSquares/h - (aver/h)*(aver/h)))* (double)h/(h - 1);
				  	
		  		if ((aver/h) != 0) e= 3.92 * (stdev/ ((Math.sqrt(h)* (aver/h))) );
		  		//e = ((Z/Math.sqrt(h))* (Math.sqrt(vrn-(aver*aver))/aver));
		  		//System.out.println("l'errore in F2 vale " + e);
		  		}
		  		i=0;
	}//fine if j
}//fine while
	
	/* STAMPA PER ERRORE
	System.out.println("COST Exited after nRUn: " +nRun);
    */
	
	//restituisco l'intervallo di confidenza di tutte le MC run
	//System.out.println( "il lower bound inter conf"+ MeanConfidenceInterval(lconf,averTot/nRun,vrn/nRun));

	averTot = averTot/nRun;
	vrnTot = (vrnTot/nRun - Math.pow(averTot, 2))*((double)nRun/(nRun-1));
	
	//System.out.println("average "+averTot);
	//System.out.println("variance  "+vrnTot);
	/*STAMPA PER ERRORE
	if(nRun==1000) System.out.println("\n CONFIDENCE COST "+MeanConfidenceInterval(lconf,averTot,vrnTot/nRun));
	*/
	
	if (nRun > 0) return MeanConfidenceInterval(lconf,averTot,vrnTot/nRun); // PASSO IL QUADRATO DELLO STANDARD ERROR
	   else return 0;
}

//Monte Carlo simulazione per il constraint che hanno parametri incerti
//sono: quello della Failure Intensity, quello sul througput del debugger (sarebbe Eq. 2 in Figura) e quello sul numero 
//dei fault detected (sarebbe Eq. 4 in Figura)

public double MonteCarloCR(double [][] C, double [] Y) {
	
	int i=0;
	int j=0;
	//variabile per memorizzare il numero di run di monte carlo in totale che eseguo
	int nRun = 0;
	//finestra di campioni che consideriamo per andare a controllare l'errore
	int h = 10;
	int k, d, t;
	double apmdk, thrg, pval, stdev;
	double aver=0.0D, sumSquares=0.0D;
	double e=0.7, fail=0.0D, app6, NTotk = 0.0, app8;
	
	//fulfillment flag binary that indicate the satisfaction of the requirements
	//boolean fReq = true;
	//tre flag
	int fReqF = 1;
	int fReqNTotk = 1;
	int fReqNTotk1 = 1;	
	int fReqDT = 1;
	
	//parametro per la somma dei flag totali
	int flagTot = 0;
	//nel caso esponenziale abbiamo solo un parametro, il beta
	double [] pSRGM = new double [K]; 
	
	//Z dobbiamo inizializzarla alla prima parte della formula dell'errore
	//1.96 e' legato al 95%
	double Z = 1.96;	
	double Tk = 0;
	int checkpval = 0;
	
	boolean check;
	double [] mu = new double [K];


	//caso di exponential
	double [] beta= new double [K];

	//qua memorizzo il risultato di ogni run (0 o 1 a seconda se sono o meno rispettati i constraint)
	List<Integer> L = new ArrayList <Integer>();
	//qua memorizzo il risultato del p-value che calcolo in funzione, di volta in volta, di tutte le MC run di monte carlo che colleziono
	List<Double> P = new ArrayList <Double>();
	
	//attenzione qua la la soglia di tolleranza e' diversa, penso non sia lo stesso valore dato nel caso
	//di MonteCarlo per stimare f[0] coie' il numero di bug corretti
	//campioniamo h valori alla volta, al h-simo andiamo a controllare l'errore
	
	//while((e > 0.005) && (nRun < 30)){
	
	while((e > 0.00001) && (nRun<1500)){
	//	System.out.println("e "+e);
	if (i < h){
		
		nRun = nRun + 1;
		//System.out.println("INIZIO calcolo");
		j=j+1;
		fail = 0;
		fReqF = 1;
		fReqNTotk = 1;
		fReqNTotk1 = 1;
		fReqDT = 1;
					
		for (k = 0; k < K; k++){	
				
				Randoms r = new Randoms();
				NTotk = 0;
				//andiamo a campionare le k funzionalita' e calcoliamo la reliability del sistema
								
				//caso di EXP usiamo uniforme
			    //passiamo alla funzione che campiona nextUniform il range di ogni parametro
			    pSRGM [k] = r.nextUniform(range[k][0], range[k][1]);
			    beta[k]= pSRGM[k];
			    
			    apmdk = mdk(k,beta[k],Y[k]);
			    
			    /*System.out.println("\n mdt nel vinc vale "+apmdk);
			    System.out.println("omega[k] vale "+ omega[k]);
			    System.out.println("Y[k] vale "+ Y[k]);
			    System.out.println("peso vale "+ weight[k]);*/
			    //System.out.println("beta vale "+ beta[k]);
			    //calcoliamo la failure intensity del sistema facendo la somma pesata della failure dei singoli servizi
				
			    //fail =0.0009;
			    fail= (weight[k]* (beta[k]*(omega[k]- apmdk))) + fail; 
			    //fail= (weight[k]* (beta[k]*(omega[k]- apmdk))) ;
			    
			    //System.out.println("fail vale "+fail);
			
			    //System.out.println("\n mdt nel vinc vale "+apmdk);
				//campioniamo un valore del numero medio di ore di lavorazione di un bug con distrib esponenziale				
				ExponentialRV Ex = new ExponentialRV(delta[k]);
				app6=Ex.nextDouble();
					//System.out.println("app6 vale "+ app6);
					//app6=10;
				
				//vedere in caso chiudere
		
				if (app6 < 0.5){
				   app6 = 1;
					}
	
			//inizio parte dei debugger
			/*	
				
	     		if (Y[k] > 0.0 && Y[k] < (B- 0.1)) 
		     		{
					   Tk = (-1/(alpha*kappa))* Math.log((Math.pow((B/Y[k]),kappa) - 1.0)/A);        	      	
				    }
					    else {if (Y[k] >= (B - 0.001)) Tk = (-1/(alpha*kappa))* Math.log((Math.pow((B/(B - 0.001)),kappa) - 1.0)/A);  
				    }
		
	
	     		for (d = 0; d < K; d++){
			
	     				if (C[d][k] >  Tk/mTh[d][k]) {							
	     						
	     					//caso in cui non e' rispettato
	     					fReqNTotk = 0;	
	     					//System.out.println("C[d][k] " +C[d][k] );
	     				//	System.out.println("C[d][k]  " +C[d][k] );
	     					//System.out.println("Tk * m vale" + mTh[d][k] * Tk);	
	     			//		System.out.println("tk vale" + Tk);
	     				}//end if flag thrgm
			
	     		}//end for sui debugger
			
	     	double somORe = 0;
	
	     	for (d = 0; d < K; d++){
	
	     		somORe = somORe + C[d][k];
		
	     		}
	
	     	if (somORe < (app6 * apmdk)) {		
	     		
	     		//caso in cui non e' rispettato
	     		fReqNTotk1 = 0;	
	     		//System.out.println("C[d][k] " +C[d][k] );
			 //   System.out.println("apmdk vale" +apmdk);
	     	}//end if flag thrg
		
	     	*/
	  
	//FINE CHIUSURA VINCOLO
	  
			    //qua considero il vincolo sul numero di fault detected (Equazione 3)		
						
				if ((apmdk - omega[k]) > 0.0){
							
						//caso in cui non e' rispettato
						fReqDT = 0;
		//				System.out.println(" yyyyy999 ");	
						}//end if flag fReqDT						
																							
	     }//end for sulle K
			
		
	    //qua vado a verificare il vincolo 1 sia rispettato
		
	   //verifichiamo che il vincolo sul Failure Rate sia rispettato
			
		if (fail > F) {
				
				//caso in cui non e' rispettato
			    
				fReqF = 0;
										
			}//end if flag fReqF
			
			//qua aggiorno il flagTot che mi dice se sono rispettati i tre requisiti
		
			flagTot = flagTot + (fReqF *  fReqNTotk *  fReqNTotk1 * fReqDT);
			  
			//memorizziamo zero se almeno un vincolo non viene rispettato, uno altrimenti
			
		  if ((fReqF *  fReqNTotk *  fReqNTotk1 * fReqDT) == 0) {
			  
			  L.add(0);					
		  			}
		  else {
				L.add(1);				
			   }
			    		
		  pval = 0;
		  //vado qua a calcolare il p-value in funzione di tutte le MC run che ho memorizzato
		  for (t = 0; t < L.size(); t++){
			  
	        	int frun = L.get(t);
	        	
	        	pval = pval + frun;
	        	//System.out.println("pval dentro value: "+pval);
	           	}
		    
		    //calcolo p-value
		    //System.out.println("pval proma value: "+pval);	
		    pval = pval/nRun;
		//      System.out.println("pval dopo value: "+pval);	
		    //  System.out.println("nrun value: "+nRun);
		    //memorizzo il p-value nella lista P  
		    
		    P.add(pval);
		   	
		   // System.out.println("List pval value: "+P);	
		   // System.out.println("il fail value: "+fail);	
		    fail = 0.0D;
			i=i+1;
		}//fine if i
	if (j >= h){
			//calcoliamo la media e la somma dei quadrati della finestra di h p-value
		  		aver=0.0D;
		  		sumSquares=0.0D;		  		
		  		//System.out.println("entro qua ");
		  		for (t = P.size()- h; t < P.size(); t++){
		  			
	        		        	double aDouble = P.get(t);
	        		        	//System.out.println("aDouble vale" + aDouble);
	        		        	aver = aver + aDouble; 
	        		        	sumSquares = sumSquares + (aDouble * aDouble);  	        		        	
	        		        	}
               
		  		if (aver != 0) {
		  			//stdev = Math.sqrt((sumSquares- aver*aver/h) / (h - 1));
		  			stdev = Math.sqrt((sumSquares/h - (aver/h)*(aver/h)))* (double)h/(h - 1);
		  			
		  			if ((aver/h) != 0)	e= 3.92 * (stdev/ ((Math.sqrt(h)* (aver/h))) );
		  		//e = ((Z/Math.sqrt(h))* (Math.sqrt(stdev-(aver*aver))/aver));
		  	    //if ((Double.isNaN(e)) || (e == Infinty)) { e=app8;}
		  		//System.out.println("l'errore di k campioni iiiii vale " + e);
		  		  }
		  		i=0;
	}//fine if j
		  	
}//fine while
		
	//ritorniamo la media dei vincoli rispettati (sarebbe il flag di robustezza)
	//calcolo la robustezza e lo restituisco
	 if (nRun > 0) return  flagTot/nRun;
	 else return 0;	
  
     
}

//si potrebbe usare la org.apache.commons.math3.stat.descriptive libreria per avere il calcolo dell'intervallo generale 
//per qualsiasi probability
/** Return the upper bound or the lower bound of the mean confidence interval  
 * and 95% approximate confidence interval*/
// nc is the number of samples
// me is the sample mean
// var is the sample variance


/**
 * Constructor.
 * Creates a default instance of the Water problem.
 * @param solutionType The solution type must "Real" or "BinaryReal".
 */

 public ReliabilityProblem97BaseIN(String solutionType) {

	//number of variables of our opt model
	//con debugger
	//numberOfVariables_  = (m * K) + K;
	 
	//senza debugger
	numberOfVariables_  =  K;
	 
	//number of obj functions
    numberOfObjectives_ = 3 ;
    
    //number of constraints
    //con debugegr
    //numberOfConstraints_= K + 2;
    //senza debugger
     numberOfConstraints_=  2;
    
    problemName_        = "ReliabilityProblem";

    lowerLimit_ = new double[numberOfVariables_];
    upperLimit_ = new double[numberOfVariables_];   
    
    //qua settiamo le variabili 
    //le prime m*K variabili menorizzano le Cdk (quindi da 0 a (m*K-1)
    //le restanti sono le Yk
     
 	//double app1=  (-1/(alpha*kappa))* Math.log((Math.pow((1),kappa) - 1.0)/A);   
 	
 	 //System.out.println("app1 vale " + app1); 
 	
    //senza debugger
    
   for (int var = 0; var < numberOfVariables_; var++){
        lowerLimit_[var] = 0.0;
        upperLimit_[var] = B;
      } //for
    
   
    
    /*
    
    //con debugger
    
    double Nmax = 100000;
    //settiamo le Cdk
    for (int var = 0; var < m*K; var++){
      lowerLimit_[var] = 0.0;
      upperLimit_[var] = Nmax;
    } //for
    
    
    //100000000 
    
    for (int var = (m*K - 1); var < numberOfVariables_; var++){
        lowerLimit_[var] = 0.0;
        upperLimit_[var] = B;
      } //for
    
   //fine con debugger 
   */
    
  //   System.out.println("Rmin vale "+ Rmin);
    
    
   
    //CODICE PER CREARE LE ISTANZE PERTURBATE 
  
    
   /*  PrintStream p;
    String s = "istanze.txt";   
    
    try {
        p = new PrintStream(s);
    }
    catch(FileNotFoundException e){
       throw new RuntimeException("errore accesso file");
    }
    
    CreazioneIstanze(p);
     p.close();
  */
              
    //questa parte indica che le soluzioni sono insiemi di variabili reali o reali e binarie
    if (solutionType.compareTo("BinaryReal") == 0)
      solutionType_ = new BinaryRealSolutionType(this) ;
    else if (solutionType.compareTo("Real") == 0)
    	solutionType_ = new RealSolutionType(this) ;
    else {
    	System.out.println("Error: solution type " + solutionType + " invalid") ;
    	System.exit(-1) ;
    } 
    
   } 
  
  /** 
  * Evaluates a solution 
  * @param solution The solution to evaluate
   * @throws JMException 
  */
  public void evaluate(Solution solution) throws JMException {
	  
    Variable[] variable  = solution.getDecisionVariables();
    
    double [] f = new double[numberOfObjectives_];
    
    double max, app;
    int d, k;
    //here we can define the variables of the opt problem
    //C are the Cdk  
    //con o senza debugger lasciare
    double [][] C = new double [m][K];
    //i number of decision variable;
    int i = 0;
  
    /*
    //togliere in caso di senza deb
    for (d = 0; d < m; d++){
      	for (k = 0; k < K; k++){
      		
      		//if (weight[k] == 0) C[d][k] = 0;
      		
      		C[d][k] = variable[i].getValue();  
      		
    	 	
    	 	if (C[d][k] > 0) x[d][k] = 1;
    	 	else x[d][k] = 0;
    	 	
    		i = i + 1;
    	}
    }
    */
    
    
    //the Y are the testing Variables
   double [] Y = new double [K];
    
    for (k = 0; k < K; k++){
       	Y[k] = variable[i].getValue();  
    	i = i+1;
    	
    }
          
    //funzione obiettivo sul numero di bug corretti
    //true in caso di lower bound
    
    // According to Max(f(x)) = -Min(f(-x)), 
    // they must be multiplied by -1. Consequently, the obtained solutions must
    // be also multiplied by -1
        
    //ATTENZIONE VEDERE IL POLONI PROBLEM
    
    f[0] = - MonteCarloF1(true, C, Y);
    
    // f[0] = -9;
    //second objective function on the time TEF inversa del testing effort 
    //andiamo a calcolare il massimo dei tempi di testing delle funzionality
   
    max=0.0;
    app=0.0;  
    //double A = 0.8; 
    //double alpha= 0.5; 
    //double kappa= 0.05;
    
    for (k = 0; k < K; k++){
   	 
      	//calcolo il tempo della funzionalita' k
    	//t_k = (-1/alpha*k) * ln[ ((B/Y_k)^(kappa) - 1)/A]
    	if (Y[k] > 0.0 && Y[k] < (B- 0.1)) 
    	{
		
			app= (-1/(alpha*kappa))* Math.log((Math.pow((B/Y[k]),kappa) - 1.0)/A);        	      	
    	
    	}
	
   else {if (Y[k] >= (B - 0.001)) app = (-1/(alpha*kappa))* Math.log((Math.pow((B/(B - 0.001)),kappa) - 1.0)/A);  
   }
       
       /* System.out.println("CALCOLO TEMPO DI TESTING funz k");
        System.out.println("la quantita' di testing allocata e' yk " + Y[k]);
        System.out.println("il tempo di testing vale " + app);*/
        
    
     if (app > max){
			max = app;
			   }
    
       }
    
    
    //vogliamo minimizzare il massimo dei tempi
    f[1] = max;
   
    //terza funzione obiettivo che minimizza i costi
    f[2]=  MonteCarloF2(false, C, Y); 
      
    
    // According to Max(f(x)) = -Min(f(-x)), 
    // they must be multiplied by -1. Consequently, the obtained solutions must
    // be also multiplied by -1
     
     
    //System.out.println("VALUAZIONE F[0] ");
    solution.setObjective(0, - 1 * f[0]);
    solution.setObjective(1, f[1]);
    solution.setObjective(2, f[2]);
    
  } // evaluate

  
  /** 
   * Evaluates the constraint overhead of a solution 
   * @param solution The solution
   * @throws JMException 
   */  
  public void evaluateConstraints(Solution solution) throws JMException {
    Variable[] variable  = solution.getDecisionVariables();
    
    double [] constraint = new double[this.getNumberOfConstraints()];
    int k;
    int d;
    int i=0;
    double app =0, totD=0, prodD;
    
    
    
  //sarebbero le N_dk  
    
    
  
  double [][] C = new double [m][K];
  //i number of decision variable;
  
//da chiudere senza debugger  
  
  /*
   for (d = 0; d < m; d++){
      	for (k = 0; k < K; k++){
      		
      		//if (weight[k] == 0) C[d][k] = 0;      		
      		C[d][k] = variable[i].getValue(); 
    	 	
      		//setto le x che sono binarie
      		
    	 	if (C[d][k] > 0) x[d][k] = 1;
    	 	else x[d][k] = 0;
    	 	
    		i = i + 1;
    	}
    }
    */
   
   
     //the Y are the testing Variables
    double [] Y = new double [K];
    
    for (k = 0; k < K; k++){
       		Y[k] = variable[i].getValue();  
    		i = i+1;
    	
    } 
      
  //scriviamo i constraint
    i=0;
    
 
  //double NOre = 0;
  
   // da togliere senza debugger
   //(k) vincoli (constraint 6 un figura)  
   /*
   for (k = 0; k < K; k++){      	
 	  
    	app = 1;
    	
    	for (d = 0; d < m; d++){ 
  	  
    		app = app * (1 - x[d][k]); 
    		   
    	}//end for su d	 
    	
    	constraint[i] = (B * (1 - app)) -  Y[k];
    	i = i+1;
    }//end for su k
  */

   // (1) vincolo sul budget totale
   app=0;
     for (k = 0; k < K; k++){      	
    	  
    	 app = Y[k] + app;
     }	  
    
     constraint[i] = B - app;   
     
 
   //  System.out.println("il vincolo vale " + constraint[i]); 
     i=i+1;
  
     
   //(mi contano come 1) constraint con parametri incerti. Sarebbero quello della Failure Intensity, quello sul througput del debugger (sarebbe Eq. 2 in Figura) 
   //e quello sul numero dei fault detected (sarebbe Eq. 4 in Figura)
   //robustness threshold (simile a qualla del paper per performance)
   
    double robT = 0.8;
     //la media dei campioni (insieme di campioni per singole funzionalita) deve essere maggiore o uguale a robT
       
   double app8= MonteCarloCR(C,Y);
       
   constraint[i] = app8 - robT; 
 
   //System.out.println("il vicolo alla fine iiiii vale " + app8);  
   //attenzione ho oscurato il calcolo della reliability del sistema
 // constraint[i] = 7;
     
    double total = 0.0;
    int number = 0;
    
    for (i = 0; i < this.getNumberOfConstraints(); i++)
    	//per i primi K-1 vincoli devono essere = 0, i va da zero..
    {
    
   	 if (constraint[i] < 0.0){
               number++;
               total+=constraint[i];
             
                              
           }
    }  
   solution.setOverallConstraintViolation(total);    
   solution.setNumberOfViolatedConstraint(number);
   //System.out.println("number alla fine iiiii vale " + number); 
   
  } // evaluateConstraints
  
 
} // ReliabilityProblem
