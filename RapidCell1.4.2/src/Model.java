import java.lang.System;
import java.io.*;
import java.util.Scanner;
public class Model {
	/**
	 * @param args
	 */
	Thread myThread;
	public Model(Thread ModelThread1){
		myThread=ModelThread1;
	}
	private int Ncells;
	public double Tmax;//s, computation time
	public double Lx, Ly, R, Xini, Yini, IniOri;// mm, length, width, radius of domain, and ini pos. of cells
	public static double commonAdaptRate;// if the population is homogeneous
	private int medium;// 0 - liquid medium, 1 - agar
	private int gradientShape;
	private double gradMaxVal, gradMinVal, gradRate, gradCenterX, gradCenterY, gradSigma;
	private double stepGradientValue, stepAdditionTime, stepRemovalTime;// parameters for gradient #9
	public Cell2D[] arrCells;
	private int boundaryCondition;//0 if reflective, 1 if periodic
	private int tumbleModel;// 0 - isotropic tumble, 1 - anisotropic tumble, angle depends on N(CW) motors
	private int runModel;// 0 - motor voting model, 1 - distortion model; 
	public static double dt;//s, integration dt 
	private double outputDt;//s, the output dt
	private int Niter=1;// inactive, used for multiple runs
	private long seed;
	private MersenneTwister RG2;
	private BufferedWriter outputStream = null, outputStreamAvg=null;
	private String fave, ffull;
	private boolean flagOutX,flagOutY,flagOutOri,flagOutCheY,//flags indicate what data goes to output files
					flagOutCheA,flagOutMeth,flagOutMb,flagOutSwFreq, 
					flagOutNmot, flagOutRun, flagOutAttr, flagOutPon;
	public long exectime = 0;
	public double SField(double x, double y, double t, int gradientShape){
		double Csignal;
		switch(gradientShape){
		case 0: return 0;
		case 1: // Constant-activity radial, with Minimum at (Lx/2,Ly/2), from 0 to 1000 Kd
			Csignal=0.999/R*(Network.Ka_on-Network.Ka_off)/Network.Kd;
			return Network.Kd*Csignal*Math.sqrt((x-Lx/2.0)*(x-Lx/2.0)+(y-Ly/2.0)*(y-Ly/2.0))/((Network.Ka_on-Network.Ka_off)/Network.Kd-Csignal*Math.sqrt((x-Lx/2.0)*(x-Lx/2.0)+(y-Ly/2.0)*(y-Ly/2.0)));
		case 2: //Constant-activity linear, along X axis, from 0 to 1000 Kd
			Csignal=0.999/Lx*(Network.Ka_on-Network.Ka_off)/Network.Kd;
			return Network.Kd*Csignal*x/((Network.Ka_on-Network.Ka_off)/Network.Kd-Csignal*x);
		case 3: //Constant-activity time ramp, from 0 to 10,000 Kd
			Csignal=0.9999/Tmax*(Network.Ka_on-Network.Ka_off)/Network.Kd;
			return Network.Kd*Csignal*(t-0)/((Network.Ka_on-Network.Ka_off)/Network.Kd-Csignal*(t-0));
		case 4: //Linear, along X
			return (gradMaxVal-gradMinVal)*x/Lx;
		case 5: //Gaussian radial, along R
			double r=Math.sqrt((x-gradCenterX)*(x-gradCenterX)+(y-gradCenterY)*(y-gradCenterY));
			return gradMaxVal*Math.exp(-r*r/(2.0*gradSigma*gradSigma));
		case 6: //Gaussian, along X
			return gradMaxVal*Math.exp(-(x-gradCenterX)*(x-gradCenterX)/(2.0*gradSigma*gradSigma));
		case 7: // exponential, in time
			//if(t<200) return 0.31*Network.Kd;// simulation from Vladimirov et al., 2008
			//else return 0.31*Network.Kd*Math.exp(5e-3*(t-200));// simulation from Vladimirov et al., 2008
			return gradRate*Math.exp(t);
		case 8: // exponential, in space
			return gradRate*Math.exp(x);
		case 9:
			if((t>=stepAdditionTime)&&(t<stepRemovalTime)) {return stepGradientValue;} // add and remove stepGradientValue [mM] of attractant
			else return 0;
		default: System.out.println("Error: Gradient type is unknown!\n"); return 0; 
		}
	}

	public void runModel() throws InterruptedException {
		exectime=java.lang.System.currentTimeMillis();
        try {
            outputStream = new BufferedWriter(new FileWriter(ffull));
            outputStreamAvg = new BufferedWriter(new FileWriter(fave));
        } catch (Exception e) {
        	System.out.println("Can not open out-files!" + e.getMessage());
          };
        	//Network.setDissConstants(1.7e-3, 12e-3, 100, 1e+6); // Aspartate binding constants
        	//Network.setDissConstants(0.02, 0.5, 100, 1e+6);  //  MeAsp  binding constants
		long nsteps=Math.round(Tmax/dt);//number of time steps
		long tsliceInterval=Math.round(outputDt/dt);//number of steps between output points
		double xprev, yprev, xnext, ynext, Snext;
// Iterations loop starts:
	for(int iter=0;iter<Niter;iter++){
		createCells(Ncells);
// Run the cells from t=0 till t=Tmax:
		for(int i=0;i<nsteps;i++){
			//if(i%100==0) {System.out.print(progressSymbols[ progressCycle ] + "\b" ); progressCycle = ( progressCycle + 1 ) & 3;	}
			if(!myThread.isInterrupted()){// listen if the STOP button is pressed
			for(int j=0;j<Ncells;j++){
				xprev=arrCells[j].getPositionX(); 
				yprev=arrCells[j].getPositionY();
				xnext=xprev+arrCells[j].cellVelocity*dt*Math.cos(arrCells[j].getOrientation()); 
				ynext=yprev+arrCells[j].cellVelocity*dt*Math.sin(arrCells[j].getOrientation());
//				Test boundary conditions:
				boolean cellReachedWall=reachedBoundary(arrCells[j], xprev, xnext, yprev, ynext, boundaryCondition);
//				Normal behavior (inside the domain):
				if(!cellReachedWall) 
				{
					Snext=SField(xnext,ynext,i*dt,gradientShape); 
					arrCells[j].chemotaxisNetwork.updateMWCmodel(Snext);
					// RUN state:
					if(arrCells[j].RunTumbleState==Cell.Run) { // If the cell runs
						arrCells[j].addRunInterval();
						arrCells[j].addRotDiffusion();
						if(medium==1){
							tryStop(arrCells[j]);// Simulation of agar traps, exponentially distributed through time
						}
						if(arrCells[j].RunTumbleState!=Cell.Stop) // works for both medium=0,1
						{
							boolean startTumble;
							if(runModel==0) startTumble=tryTumbleVoting(arrCells[j]); //voting model
							else	startTumble=tryTumbleDistortionModel(arrCells[j]);//cw-distortion model
							if(startTumble) {
								double NMotorsCW=arrCells[j].getNMotorsCW();
								double angle;
								if(tumbleModel==0) angle=arrCells[j].tumble();// use System.out.print() to get the angle
								else angle=arrCells[j].tumbleCWdependent(NMotorsCW);
							};
						}
					}
					//TUMBLE state:
					else if(arrCells[j].RunTumbleState==Cell.Tumble) { //If the cell tumbles
						arrCells[j].addTumbleInterval();
						//voting model:
						if (runModel==0){ 
						xnext=xprev; ynext=yprev;
						tryRunVoting(arrCells[j]);
						}
						//cw-distortion model:
						else tryRunDistortionModel(arrCells[j]);
					}
					// STOP state (in agar only):
					else if(arrCells[j].RunTumbleState==Cell.Stop) {
						//voting model:
						if (runModel==0){ 
							xnext=xprev; ynext=yprev;
							tryTumbleVoting(arrCells[j]);
							}
							//cw-distortion model:
							else tryTumbleDistortionModel(arrCells[j]);
					}
					updateMotorStates(arrCells[j]); //update CW/CCW states of the motors.
					arrCells[j].setPosition(xnext,ynext);// update cell position.	
				}
			} // end of cycle for(int j=0;j<Ncells;j++)

			//Record the cells states:
			if(((i%tsliceInterval)==0)||(i==nsteps-1))
			{
			double sumPositionX=0, sumPositionY=0, sumCheYp=0, sumOri=0, sumCheA=0, sumMeth=0, sumMb=0, 
			sumSwFreq=0, sumRun=0, sumAttr=0, sumNmot=0, sumPon=0;
			try {
				outputStream.write(String.format("%3.2f\t", i*dt));// current time point
				for(int j=0;j<Ncells;j++){		
					if(flagOutX) outputStream.write(String.format("%2.4f\t", arrCells[j].getPositionX()));
					if(flagOutY) outputStream.write(String.format("%2.4f\t", arrCells[j].getPositionY()));
					if(flagOutOri) outputStream.write(String.format("%2.4f\t", arrCells[j].getOrientation()));
					if(flagOutCheA) outputStream.write(String.format("%1.4f\t", arrCells[j].chemotaxisNetwork.cheAp));
					if(flagOutCheY) outputStream.write(String.format("%1.4f\t", arrCells[j].chemotaxisNetwork.cheYp));
					//outputStream.write(arrCells[j].chemotaxisNetwork.adaptRate+"\t");
					if(flagOutMeth) outputStream.write(String.format("%1.4f\t", arrCells[j].chemotaxisNetwork.meth));
					//outputStream.write(arrCells[j].chemotaxisNetwork.CheR+"\t");
					//outputStream.write(arrCells[j].chemotaxisNetwork.CheBP+"\t");
					if(flagOutMb) outputStream.write(String.format("%1.4f\t", arrCells[j].flagellarMotors[0].getMotorBias()));	
					if(flagOutSwFreq) outputStream.write(String.format("%1.4f\t", arrCells[j].flagellarMotors[0].getCCWtoCWswitchFreq()));
					if(flagOutNmot)	outputStream.write(String.format("%1.0f\t", arrCells[j].getNMotorsCW()));
					if(flagOutRun)  outputStream.write(String.format("%1.0f\t", arrCells[j].RunTumbleState));
					double xattr=arrCells[j].getPositionX();
					double yattr=arrCells[j].getPositionY();
					if(flagOutAttr){
						outputStream.write(String.format("%1.6f\t", SField(xattr,yattr,i*dt,gradientShape)));
					}
					if(flagOutPon)  outputStream.write(String.format("%1.4f\t", arrCells[j].chemotaxisNetwork.P_on));
					// Compute the average values:
					sumPositionX+=arrCells[j].getPositionX();
					sumPositionY+=arrCells[j].getPositionY();
					sumCheYp+=arrCells[j].chemotaxisNetwork.cheYp;
					sumOri+=arrCells[j].getOrientation();
					sumCheA+=arrCells[j].chemotaxisNetwork.cheAp;
					sumMeth+=arrCells[j].chemotaxisNetwork.meth;
					sumMb+=arrCells[j].flagellarMotors[0].getMotorBias();
					sumSwFreq+=arrCells[j].flagellarMotors[0].getCCWtoCWswitchFreq();
					sumRun+=arrCells[j].RunTumbleState;
					sumAttr+=SField(xattr,yattr,i*dt,gradientShape);
					sumNmot+=arrCells[j].getNMotorsCW();
					sumPon+=arrCells[j].chemotaxisNetwork.P_on;
					}
			outputStream.write("\n");
			outputStreamAvg.write(String.format("%3.2f\t",i*dt));
			if(flagOutX) outputStreamAvg.write(String.format("%2.4f\t",sumPositionX/Ncells));
			if(flagOutY) outputStreamAvg.write(String.format("%2.4f\t",sumPositionY/Ncells));
			if(flagOutOri) outputStreamAvg.write(String.format("%2.4f\t",sumOri/Ncells));
			if(flagOutCheA) outputStreamAvg.write(String.format("%2.4f\t",sumCheA/Ncells));
			if(flagOutCheY) outputStreamAvg.write(String.format("%2.4f\t",sumCheYp/Ncells));
			if(flagOutMeth) outputStreamAvg.write(String.format("%2.4f\t",sumMeth/Ncells));
			if(flagOutMb) outputStreamAvg.write(String.format("%2.4f\t",sumMb/Ncells));
			if(flagOutSwFreq) outputStreamAvg.write(String.format("%2.4f\t",sumSwFreq/Ncells));
			if(flagOutNmot) outputStreamAvg.write(String.format("%2.4f\t",sumNmot/Ncells));
			if(flagOutRun) outputStreamAvg.write(String.format("%2.4f\t",sumRun/Ncells));
			if(flagOutAttr) outputStreamAvg.write(String.format("%2.4f\t",sumAttr/Ncells));
			if(flagOutPon) outputStreamAvg.write(String.format("%1.4f\t",sumPon/Ncells));
			outputStreamAvg.write("\n");
			}
			catch (IOException e) {
				System.out.println(e.getMessage());
			};
			}
			}// end of if(!myThread.isInterrupted()){
		}// END of i=(0..nsteps) time loop
	}//END of "iter" loop 
	
	exectime=java.lang.System.currentTimeMillis()-exectime;
    try {
		outputStream.close();
		outputStreamAvg.close();
	} catch (IOException e) {
		System.out.println(e.getMessage());
	}   
	}// end of main()
////////////////////////////////////////////////////////
///////////////////FUNCTIONS////////////////////////////	
////////////////////////////////////////////////////////
	
public void createCells(int newNcells){
		Ncells=newNcells;
		arrCells= new Cell2D[Ncells];
		for(int j=0;j<Ncells;j++){
			double x = Xini, y = Yini;

			// Random initial position
			if (x == -1) x = Lx * RG2.nextDouble();
			if (y == -1) y = Ly * RG2.nextDouble();

			if(IniOri==1){
				arrCells[j]=new Cell2D(x,y);// random ori
			}
			else{
				arrCells[j]=new Cell2D(x,y,IniOri);
			}
		}

	}

/////////////////////////////////////////////////////////
public boolean tryTumbleDistortionModel(Cell cell_j){
	boolean startTumble=false;
	int Nmotors=Cell.Nmotors;
	double NrotatingCW=0, sumCW=0;
	for(int j=0;j<Nmotors;j++){
		if((cell_j.flagellarMotors[j].getCWintervalDuration()>0))
		{	NrotatingCW++;	}
		if(cell_j.flagellarMotors[j].getCWintervalDuration()<=0.15)
		{
			sumCW=sumCW+(cell_j.flagellarMotors[j].getCWintervalDuration()/0.15);
		}
		else if(cell_j.flagellarMotors[j].getCWintervalDuration()>0.15)
		{
			sumCW=sumCW+Math.exp(-15.0*(cell_j.flagellarMotors[j].getCWintervalDuration()-0.15));
		}
		
	}
	cell_j.cellVelocity=Cell.maxCellVelocity*Math.max(0,1.0-sumCW);
	//System.out.print(cell_j.cellVelocity+"\n");
	if(sumCW>=1.0) {// 
		startTumble=true;
		cell_j.RunTumbleState=Cell.Tumble;
		//System.out.print(cell_j.getRunIntervalDuration()+"\t");
		cell_j.resetRunIntervalDuration();
	}
	return startTumble;
}

/////////////////////////////////////////////////////////
	public boolean tryTumbleVoting(Cell cell_j){
		boolean startTumble=false;
		int Nmotors=Cell.Nmotors;
		double NrotatingCW=0;
		for(int j=0;j<Nmotors;j++){
			if((cell_j.flagellarMotors[j].getCWintervalDuration()>0))
			{	NrotatingCW++;	}
		}
		if((NrotatingCW/Nmotors)>0.5) {// start a new tumble
			startTumble=true;
			cell_j.RunTumbleState=Cell.Tumble;
			//System.out.print(cell_j.getRunIntervalDuration()+"\n");
			cell_j.resetRunIntervalDuration();
		}
		return startTumble;
	}
/////////////////////////////////////////////////////////
	public boolean tryRunDistortionModel(Cell cell_j){
		boolean startRun=false;
		int Nmotors=Cell.Nmotors;
		double NrotatingCCW=0, sumCW=0;
		for(int j=0;j<Nmotors;j++){
			if((cell_j.flagellarMotors[j].getCCWintervalDuration()>0))
			{	NrotatingCCW++; 	}
			if(cell_j.flagellarMotors[j].getCWintervalDuration()<=0.15)
			{
				sumCW=sumCW+(cell_j.flagellarMotors[j].getCWintervalDuration()/0.15);
			}
			else if(cell_j.flagellarMotors[j].getCWintervalDuration()>0.15)
			{
				sumCW=sumCW+Math.exp(-20.0*(cell_j.flagellarMotors[j].getCWintervalDuration()-0.15));
			}
		}
		cell_j.cellVelocity=Cell.maxCellVelocity*Math.max(0,1.0-sumCW);
		if(sumCW<1.0) {// Stop tumbling, run:
			cell_j.RunTumbleState=Cell.Run; 
			startRun=true;
			//System.out.print(arrCells[j].getTumbleIntervalDuration()+"\n");
			cell_j.resetTumbleIntervalDuration();
			}
		return startRun;
	}

/////////////////////////////////////////////////////////
	public boolean tryRunVoting(Cell cell_j){
		boolean startRun=false;
		int Nmotors=Cell.Nmotors;
		double NrotatingCCW=0;
		for(int j=0;j<Nmotors;j++){
			if((cell_j.flagellarMotors[j].getCCWintervalDuration()>0))
			{	NrotatingCCW++; 	}
		}
		if(((NrotatingCCW/Nmotors)>0.5)) {// Stop tumbling, run:
			cell_j.RunTumbleState=Cell.Run; 
			startRun=true;
			//System.out.print(cell_j.getTumbleIntervalDuration()+"\n");
			cell_j.resetTumbleIntervalDuration();
			}
		return startRun;
	}
/////////////////////////////////////////////////////////
	public void tryStop(Cell cell_j){
		double rand01=RG2.nextDouble();
		double velocity=20e-3;
		double meanFreePath=40e-3; //mm, mean free path between collisions 
		double stopFreq=velocity/meanFreePath;
		if(rand01<=dt*stopFreq) {cell_j.stop();}
	}
/////////////////////////////////////////////////////////
	public void updateMotorStates(Cell cell_j){
		double rand01,rand02;
		// Multiple motors
		int Nmotors=Cell.Nmotors;
		for(int j=0;j<Nmotors;j++){
			rand01=RG2.nextDouble();
// Try switch CCW->CW:
			if(cell_j.flagellarMotors[j].RotaryState==1){
				cell_j.flagellarMotors[j].addDtToCCWinterval();
				// Try to switch CCW->CW:
				if((rand01<=cell_j.flagellarMotors[j].getCCWtoCWswitchFreq()*dt))
				{ // switch to CW-rotation
					cell_j.flagellarMotors[j].RotaryState=0;
					//System.out.print(cell_j.flagellarMotors[j].getCCWintervalDuration()+"\n");
					cell_j.flagellarMotors[j].resetCCWintervalDuration();
				}
			}
			rand02=RG2.nextDouble();
// Try switch CW->CCW:
		if(cell_j.flagellarMotors[j].RotaryState==0){
			cell_j.flagellarMotors[j].addDtToCWinterval();
			if((rand02<=cell_j.flagellarMotors[j].getCWtoCCWswitchFreq()*dt)){
				cell_j.flagellarMotors[j].RotaryState=1;
				//System.out.print(cell_j.flagellarMotors[j].getCWintervalDuration()+"\n");
				cell_j.flagellarMotors[j].resetCWintervalDuration();
			}
		}
	}
	}
/////////////////////////////////////////////////////////
	private boolean reachedBoundary(Cell cell_j, double xprev, double xnext, double yprev, double ynext, int boundaryCondition0){
		boolean reachedWall=false;
		//		Plain	walls:
		if(boundaryCondition0==0){// reflective boundary
			if(xnext<0) { 
				cell_j.reverseOrientation(0); 
				cell_j.tumble(); 
				reachedWall=true;
				}
			else if(xnext>Lx) { 
				cell_j.reverseOrientation(2); 
				cell_j.tumble(); 
				reachedWall=true;
				}
			else if(ynext<0) { 
				cell_j.reverseOrientation(3); 
				cell_j.tumble(); 
				reachedWall=true;
				}
			else if(ynext>Ly) { 
				cell_j.reverseOrientation(1); 
				cell_j.tumble(); 
				reachedWall=true;
				}
		}
		else{// periodic boundary
			if(xnext<0) { 			
				cell_j.setPosition(Lx+xnext, ynext);
				reachedWall=true;
			}
			else if(xnext>Lx) { 
				cell_j.setPosition(xnext-Lx, ynext);
				reachedWall=true;
			}
			else if(ynext<0) { 
				cell_j.setPosition(xnext, Ly+ynext);
				reachedWall=true;
			}
			else if(ynext>Ly) { 
				cell_j.setPosition(xnext, ynext-Ly);
				reachedWall=true;
			}
		}
		return reachedWall;

	}
/////////////////////////////////////////////////////////
	// SETTERS & GETTERS //
	public void setNcells(int Ncells1){	Ncells=Ncells1;};
	public void setTmax(double Tmax1){	Tmax=Tmax1;};
	public void setLxLy(double Lx1, double Ly1){	Lx=Lx1; Ly=Ly1; R=Lx1/2.0;};
	public void setInitialCellDistribution(double Xini0, double Yini0, int OriIndex){
		Xini=Xini0;
		Yini=Yini0;
		switch (OriIndex){//select the initial orientation of cells
		case 0: IniOri=0; 			break;
		case 1: IniOri=Math.PI/2.0; break;
		case 2: IniOri=Math.PI; 	break;
		case 3: IniOri=3*Math.PI/4;	break;
		case 4: IniOri=1;			break; // for random ori, works like flag
		default: IniOri=1; 			break;
		}
		
		};
	public void setGradient(int Grad){gradientShape=Grad;};
	public void setMedium(int Med){medium=Med;};
	public void setRunModel(int rmod){runModel=rmod;};
	public void setTumbleModel(int tmod){tumbleModel=tmod;};
	public void setAdaptRate(double arate){commonAdaptRate=arate;};
	public void setStepGradientParams(double stepConc, double t_add, double t_rem){
		stepGradientValue=stepConc;
		stepAdditionTime=t_add;
		stepRemovalTime=t_rem;
		}
	public void setGradMinMaxVal(double min1, double max1){
		gradMinVal=min1;
		gradMaxVal=max1;
	}
	public void setGradSourceXY(double x, double y){
		gradCenterX=x;
		gradCenterY=y;
	}
	public void setExpGradRate(double rate){gradRate=rate;}
	public void setGradSigma(double sigma){gradSigma=sigma;}
	public void setOutputDt(double dt0){outputDt=dt0;};// sets the output time step
	public void setOutFileNames(String fave0, String ffull0){//sets the file names for output
		fave=fave0;
		ffull=ffull0;
	}
	public void setOutFlags(boolean flagOutX0, boolean flagOutY0, boolean flagOutOri0, boolean flagOutCheY0,
			boolean flagOutCheA0, boolean flagOutMeth0, boolean flagOutMb0, boolean flagOutSwFreq0, 
			boolean flagOutNmot0, boolean flagOutRun0, boolean flagOutAttr0, boolean flagOutPon0){
		flagOutX=flagOutX0;
		flagOutY=flagOutY0;
		flagOutOri=flagOutOri0;
		flagOutCheY=flagOutCheY0;//flags indicate what data goes to output files
		flagOutCheA=flagOutCheA0;
		flagOutMeth=flagOutMeth0;
		flagOutMb=flagOutMb0;
		flagOutSwFreq=flagOutSwFreq0;
		flagOutNmot=flagOutNmot0;
		flagOutRun=flagOutRun0;
		flagOutAttr=flagOutAttr0;
		flagOutPon=flagOutPon0;
	}
	public void setDt(double dt0){dt=dt0;}
	public void setRandSeed(long s){
		seed=s;
		RG2 = new MersenneTwister(seed);
		}
	public void setBoundaryCondition(int newBoundaryCond){boundaryCondition=newBoundaryCond;}
} // end of class Model 
