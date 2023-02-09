package jmetal.core;

import java.util.ArrayList;


public class AdditionalInfo {

	private ArrayList<Double> costVector;
	private ArrayList<Double> timeVector;
	private long time;

	public long getTime() {
		return time;
	}

	public void setTime(long time) {
		this.time = time;
	}

	public AdditionalInfo(){
		costVector = new ArrayList<Double>();
		timeVector = new ArrayList<Double>();
	}
	
	public void setCostVector(ArrayList<Double> _costVector) {
		// TODO Auto-generated method stub
		costVector = _costVector;
	}

	public void setTimeVector(ArrayList<Double> _timeVector) {
		timeVector = _timeVector;
	}

	public ArrayList<Double> getCostVector() {
		return costVector;
	}

	public ArrayList<Double> getTimeVector() {
		return timeVector;
	}

}
