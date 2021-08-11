import java.util.Arrays;
import java.util.LinkedList;
public class Network {
	protected double dt=Model.dt;
	public double cheYp;
	private static double CheA_tot, CheR, CheBP, CheY_tot, CheZ, meth0;
	private static final double CheR_wt=1, CheB_wt=1;// 
	// parameters of MWC model:
	public static double Ka_off;// Tar, mM, 0.02 - MeAsp[Endres06], 1.7e-3 - Asp [Emonet05]; 0.014 - Tar, MeAsp[Sourjik, Temp. effects, unpublished]
	public static double Ka_on;// Tar, mM, 0.5 - MeAsp[Endres06], 12e-3 - Asp [Emonet05]
	public static double Kd=Math.sqrt(Ka_off*Ka_on);//
	private static double Ks_off;// Mm, 100 - Tsr, MeAsp[Endres06]
	private static double Ks_on;// mM, 1e+6 - Tsr, MeAsp[Endres06]
	private static double n_a;// 6 Tars in a cluster of 18
	private static double n_s;// 12 Tsrs in cluster of 18
	public double meth;// methylation state (0..8)
	public double cheAp;
	public double P_on;
	public double adaptRate;// adaptation rate
	public int groupID=0;// grouping variable, to simulate heterogeneous populations
	private static double adaptPrecision;// set 1.0 for a perfect adaptation (100%)
	private static double k_R=0.0182, k_B=0.0364;// effective catalytic rates
	private static double kA=5, kY=100, kZ=30, gY=0.1;// 
	private static long seed;
	private static MersenneTwister RG3;
	private LinkedList<Double> methMemory = new LinkedList<Double>();
	private static int ligandMethylationRelationIndex;
	private int steps;
	private final int stepsPerMethChange = 1000;

	public Network(Cell newCell){
		adaptRate=Model.commonAdaptRate;//default coefficient of adaptation rate, if the population is homogeneous
		meth=meth0;
		P_on=1./3.;
	}

	private double getNextExponential(double lambda) {
	    return  Math.log(1 - RG3.nextDouble()) / (-lambda);
	}

	private double ligandToMethPoly(double logS) {
		if (logS > 3.3) return 8.0;

		return 5.150327199838268
			   + 1.1960663826447897 * Math.pow(logS, 1)
			   - 0.4712982846490096 * Math.pow(logS, 2)
			   - 0.04125276995184807 * Math.pow(logS, 3)
			   + 0.08202995456148689 * Math.pow(logS, 4)
			   + 0.00012251899823295862 * Math.pow(logS, 5)
			   - 0.0038056894330643696 * Math.pow(logS, 6)
			   + 0.00015745073841346674 * Math.pow(logS, 7);
	}

	private void clampMeth() {
		if (meth < 0.0) meth = 0.0;
		else if (meth > 8.0) meth = 8.0;
	}

	private void updateMeth() {
		// update methylation level, by simple ordinary differential equation
		meth = meth + dt * adaptRate * (k_R * CheR * (1.0 - P_on) - k_B * CheBP * P_on);
	}

	public void updateMWCmodel(double S){
		double eps_meth = 0.0, sum_fa, sum_fs, F, logS = Math.log10(S);
		steps++;

		// TODO: add random gaussian noise to precise function
		switch (ligandMethylationRelationIndex) {
			case 0:  // Real methylation dynamics (maximum drift)
				updateMeth();
				break;
			case 1:  // Step function (high drift but with reduced information)
				meth = Math.round(4 * ligandToMethPoly(logS)) / 4.0;
				break;
			case 2:  // Gaussian distribution (medium drift with minimal information)
				if (steps % stepsPerMethChange == 0) meth = 4.0 + RG3.nextGaussian();
				break;
			case 3:  // M = 0.0 (minimum entropy and minimum information)
				meth = 0.0;
				break;
			case 4:  // Exponential distribution (very low entropy and information with some drift)
				if (steps % stepsPerMethChange == 0) meth = getNextExponential(0.5);
				break;
			case 5:  // Only includes two linearly increasing segments, otherwise m = 0.0 (balance of drift and entropy, leaning towards maximizing drift)
				if ((logS < -1.5) || ((logS > 1.0) && (logS < 2.5)) || (logS > 3.25)) meth = 0.0;
				else {
					if (meth == 0.0) meth = ligandToMethPoly(logS);
					updateMeth();
				}
				break;
			case 6:  // Only includes one linearly increasing segment, otherwise m = 0.0 (balance of drift and entropy, leaning towards minimizing entropy)
				if ((logS < 2.5) || (logS > 3.25)) meth = 0.0;
				else {
					if (meth == 0.0) meth = ligandToMethPoly(logS);
					updateMeth();
				}
				break;
			case 7:  // 7th degree polynomial best fit based on 100 simulations
				meth = ligandToMethPoly(logS);
				break;
			// case 8:  // Memory of previous ligand levels using S from several time steps ago
			// 	methMemory.add(S);
			// 	if (methMemory.size() > 100) S = methMemory.remove();
			// 	else S = methMemory.element();
			// 	updateMeth();
			// 	break;
			// case 9:  // Piecewise linear methylation level as a function of ligand concentration
			// 	if (logS <= -2.0) meth = 1.85;
			// 	else if (logS <= 1.0) meth = 4.0 / 3.0 * logS + 14.0 / 3.0;
			// 	else if (logS <= 2.0) meth = 0.5 * logS + 5.5;
			// 	else meth = 1.5 * logS + 3.5;
			// 	meth = meth + 0.1 * RG3.nextGaussian();
			// 	meth = Math.max(0.0, Math.min(8.0, meth));
			// 	break;
		}

		// Clamp methylation level to be between 0 and 8
		clampMeth();

		// compute free-energy offset from methylation:
		if(meth<0) eps_meth=1.0; // Endres & Wingreen, 2006, piece-wise linear
			else if(meth<2) eps_meth= 1.0-0.5*meth;
				else if (meth<4) eps_meth=-0.3*(meth-2.0);
					else if  (meth<6) eps_meth=-0.6-0.25*(meth-4.0);
						else if (meth<7) eps_meth=-1.1-0.9*(meth-6.0);
							else if (meth<8) eps_meth=-2.0-(meth-7.0);
								else eps_meth=-3.0;
		sum_fa= n_a*(eps_meth + Math.log((1+S/Ka_off)/(1+S/Ka_on)));
		sum_fs= n_s*(eps_meth + Math.log((1+S/Ks_off)/(1+S/Ks_on)));
		F=sum_fa+sum_fs;
		P_on=1.0/(1.0+Math.exp(F)); //probability that rec. cluster is ON
		cheAp=CheA_tot*P_on*kA/(P_on*kA+kY*CheY_tot); //amount of phosphorylated CheA
		double scaling=19.3610;// scales cheYp linearly so that cheYp=1 at rest (p_on=1/3)
		cheYp=adaptPrecision*scaling*CheY_tot*kY*cheAp/(kY*cheAp+kZ*CheZ+gY);// phosphorylated CheY-P, scaled to 1.0 at rest
		// double dFdc = n_a * (S / (S + Ka_off) - S / (S + Ka_on)) + n_s * (S / (S + Ks_off) - S / (S + Ks_on));
		// double drift = 1.0 * P_on * (1 - P_on) * dFdc * 0.75;
		// System.out.println("Drift = " + drift);
	}
	public static void setNreceptors(int Ntar, int Ntsr){
		n_a=Ntar;// Tars in a cluster 
		n_s=Ntsr;// Tsrs in a cluster
	}
	public void setAdaptRate(double newAdaptRate){
		adaptRate=newAdaptRate;
	}

	public static void setDissConstants(double Ka_off1, double Ka_on1, double Ks_off1, double Ks_on1){
		Ka_off=Ka_off1;
		Ka_on=Ka_on1;
		Kd=Math.sqrt(Ka_off*Ka_on);
		Ks_off=Ks_off1;
		Ks_on=Ks_on1;
	}
	public static void setAdaptPrecision(double AP){adaptPrecision=AP;}
	public static void setRelativeCheA(double A) {CheA_tot = A;}
	public static void setRelativeCheRCheB(double R, double B){	CheR=R*CheR_wt; CheBP=B*CheB_wt;}
	public static void setRelativeCheYCheZ(double Y, double Z){	CheY_tot=Y; CheZ=Z;}
	public static void setIniMethState(double m){meth0=m;}
	public static void setRates(double kR0, double kB0, double kA0, double kY0, double kZ0){
		k_R=kR0;
		k_B=kB0;
		kA=kA0;
		kY=kY0;
		kZ=kZ0;
	}
	public static void setLigandMethylationRelation(int relationIndex) {ligandMethylationRelationIndex = relationIndex;}
	public static void setRandSeed(long s){
		seed = s;
		RG3 = new MersenneTwister(seed);
	}
}
