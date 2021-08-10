public abstract class Cell {
	public static int lastID=0;
	public int CellID;
	public double RunTumbleState;
	public static final double Tumble=0;
	public static final double Run=1;
	public static final double Stop=2;
	public double RunIntervalDuration=0;
	public double TumbleIntervalDuration=0;
	protected double dt=Model.dt;
	public Network chemotaxisNetwork;
	public static int Nmotors;// 5 by default
	public Motor[] flagellarMotors;
	public static double Drotation;//rad^2/s, 0.062 [Berg, 1983]; 0.1605 [Andrews,2006, T=20C]; 0.2086 (T=30C)
	public static double maxCellVelocity;// mm/s, maximum cell velocity
	public double cellVelocity=maxCellVelocity;// mm/s, actual cell velocity
	public double MaxCWmotorsInTumble=0;// maximal number of motors involved in the current tumble
	//Normal cell:
	public Cell(){
		chemotaxisNetwork=new Network(this);
		flagellarMotors=new Motor[Nmotors];
		for(int j=0;j<Nmotors;j++){ flagellarMotors[j]=new Motor(this); }
	}
	public void addRunInterval(){
		RunIntervalDuration=RunIntervalDuration+dt;
	}
	public void addTumbleInterval(){
		TumbleIntervalDuration=TumbleIntervalDuration+dt;
	}
	public void resetRunIntervalDuration(){
		RunIntervalDuration=0;
	}
	public void resetTumbleIntervalDuration(){
		TumbleIntervalDuration=0;
	}
	public double getRunIntervalDuration(){
		return RunIntervalDuration;
	}	
	public double getTumbleIntervalDuration(){
		return TumbleIntervalDuration;
	}	
	public abstract double tumble();
	public void stop(){
		RunTumbleState=Stop;
	}
	public abstract double getPositionX();
	public abstract double getPositionY();
	public abstract void setPosition(double newPosX,double newPosY);
	public abstract void reverseOrientation(int border);
	public abstract double getOrientation();
	public abstract void setOrientation(double newOri);
	public void addRotDiffusion(){};
	public static void setRotDiffusionCoeff(double Dr){Drotation=Dr;};
	public static double getRotDiffusionCoeff(){return Drotation;};
	public void setCellVelocity(double V){cellVelocity=V;};
	public static void setMaxCellVelocity(double V){maxCellVelocity=V;};
	public void setRunTumbleState(double state){RunTumbleState=state;};
	public double getNMotorsCW(){
		double NrotatingCW=0;
		for(int j=0;j<Nmotors;j++){
			if((flagellarMotors[j].getCWintervalDuration()>0))	
				{NrotatingCW++;	}
		}	
		return NrotatingCW;
	}
	public static void setNmotors(int N){Nmotors=N;};
}
