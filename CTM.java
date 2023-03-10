/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ped;


import static java.lang.Math.pow;
import util.DoubleE;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;


public class CTM extends Link {

	/**
	 * Cell[] is the array of cells used to model this link, which means each elements in the array is a Cell
	 * This array has NOT been instantiated! You need to do that.
	 */
	private Cell[] cells;

	// TODO: create two LinkedList store the cumulative N for the first cell and last cell
	//    * How to code them in a smart way
	private LinkedList<Double> N_up;
	private LinkedList<Double> N_down;

	private double yin;
	private double yout;


	public CTM(int id, Node source, Node dest, double length, double ffspd, double capacityPerLane, int numLanes) {
		super(id, "CTM", source, dest, length, ffspd, capacityPerLane, numLanes);


		// TODO: read in length() units is ft ---> but we should use mile here!!!
		//    * Note the units !!!
		// int numCells = (int) Math.round(length / getCellLength()); // ---> This is MWL code, he read mile from txt.files
		int numCells = (int) Math.round( (length / 5280) / getCellLength()); // ---> you should change unitsg

		cells = new Cell[numCells]; // allocate the size of cells array

		/**
		 * NOTE: the first index is "0"
		 */
		for (int i = 0; i < cells.length; i++) {
			cells[i] = new Cell(this);
		}
		// Notice the index start from
		for (int i = 1; i < cells.length; i++) {
			cells[i].setPrevCell(cells[i - 1]);
			cells[i - 1].setNextCell(cells[i]);
		}

		reset();


	}




	public void reset() {

		N_up = new LinkedList<>();
		N_down = new LinkedList<>();

		for (Cell c : cells) {
			c.reset();
		}

		/**
		 * add first cell and last cell cumulative N
		 */


		for (int i = 0; i < Math.round(getLength()/5280 / (getFFSpeed() / 3600.0) / Params.dt) + 1; i++) {
			N_up.add(0.0);
		}

		System.out.println( Math.round(getLength()/5280 / (getFFSpeed() / 3600.0) / Params.dt) + 1);

		// TODO: When running Sioux-falls, file reader can read in w_b from txt file directly
		//     *: If change from Triangular-FD to Trapezoidal-FD, q_max should be lower than before?
		//     *: we may not use original q_max, when we want to simulate CTM based Trapezoidal-Fd
		double w =  getBackwardSpeed(); // ---> This is defined in Link.java this is used for My'APBWP checking

		for (int i = 0; i < Math.round(getLength()/5280 / (w / 3600.0) / Params.dt) + 1; i++) {
			N_down.add(0.0);
		}


		//System.out.println("N_down initialize testing " + getId() + " " + N_down.size());
		System.out.println("N_down initialize testing " + getId() + " " + N_down);
		System.out.println("N_down initialize testing "+ " " + getId() + " " + N_up);


		/**
		 * This is reset
		 */
		yin = 0;
		yout = 0;


	}


	// here getOccupancy is number of vehicles on this CTM.link
	public double getOccupancy(){
		double output = 0.0;
		for (Cell c: cells){
			output += c.getOccupancy();
		}
		return output;
	}


	// unit is mile/hour *  (one time-step / 3600 seconds )
	public double getCellLength()
	{
		//System.out.println ("length of each Cell is ---> " +  " " + getCellLength());
		return (getFFSpeed() * Params.dt/3600.0); // dx = uf * dt = 30 mile / hour * time-step-duration = 30 / 3600 * 6 seconds
	}

	// freeflow travel time ???
	// TODO 理解 each cell free flow travel time: dt; cell.length = number of cells; so freeflow travel time for a ctm link is ===> N * dt
	public double getFFTime()
	{
		// notice that cells = sum of Cell
		// cells.length here is the number of Cell for a CTM link
		return cells.length * Params.dt; // cells is subclass of Link.java cells.length =  return units is mile
	}



	public void step()
	{
		for(Cell c : cells)
		{
			c.step();

		}
	}


	public void update(){


		int cell_id = 0;
		for (Cell c: cells){
			System.out.println("Before: Testing each cell's number of veh for testing demand loading --->" + "current cell_id: "+ " " + cell_id + " " + " occupancy is "+ " "  + c.getOccupancy());
			c.update();
			// checking
			//System.out.println("After: Testing each cell's number of veh for testing demand loading --->" + "current cell_id: "+ " " + cell_id + " " + " occupancy is "+ " "  + c.getOccupancy());
			cell_id += 1;
		}




		// TODO: cells[0].getSendingFlow() ??? OR cells[0].getOccupancy()
		//    * getOccupancy()的问题是，如果是红灯，他会一直叠加上次的n
		//    * getSendingflow的问题也是，如果是红灯，他还是会叠加上次的的n

		//System.out.println(N_up); // ---> debug tips: print out before the bug lines

		double N_last = N_up.getLast();
		N_up.add(N_last + yin);



		N_last = N_down.getLast();
		N_down.add(N_last + yout);

		yin = 0;
		yout = 0;

		N_up.removeFirst();
		N_down.removeFirst();

		// TODO: testing_dt printout current time steps
		double current_timestep= Params.time / Params.dt;

		// TODO: 这样打印的顺序更加直观
		System.out.println(current_timestep + " " + getId() + " " + N_down);
		System.out.println(current_timestep + " " + getId() + " " + N_up);

		step();

	}


	// max flow that CTM link can serve at next time step
	public double getReceivingFlow()
	{
		return cells[0].getReceivingFlow(); // first cell
	}

	// max flow that CTM link can send to the downstream link at next time step
	public double getSendingFlow()
	{
		return cells[cells.length-1].getSendingFlow(); // last cell
	}


	public void addFlow(double y)
	{
		//yin += cells[0].addFlow(y);
		// first cell addFlow(y)
		cells[0].addFlow(y); // first cell
		logEnteringFlow(y);

		// TODO: This is for APWBP
		yin += y;
	}


	public void removeFlow(double y)
	{
		cells[cells.length-1].removeFlow(y); // last cell 每次移除的flow

		// TODO: This is for APWBP
		yout += y;
	}




	// -------------------------------------------------- My paper's Algorithm -----------------------------------------

	public LinkedList<Double> ShockWaveDetection()
	{
		LinkedList<Double> correct_x_shock = new LinkedList<Double>();

		// TODO: you can not put temp_xs as 0, since it may queue to the start end of the link
		double temp_xs;
		// TODO changing uf wb units, per/dt & per/seconds is huge difference
		double uf = ((this.getFFSpeed()/3600)*(Params.dt));
		double wb = (this.getBackwardSpeed()/3600)*(Params.dt);
		double K_jam = Params.JAM_DENSITY;

		// TODO changing length units, Yash is ft & MWL is mile
		double L = (this.getLength()/5280);
		System.out.println("uf--->"+ " " + uf + " "  + "wb--->"+ " " + wb + " " + "K_jam--->"+ " " + K_jam + " " + "L--->" + " " + L + "for this link, linkid is--->" + this.getId());

		int current_time = (int)Math.floor(Params.time/Params.dt); // 感觉向下取整比较合理 0-6 seconds 属于 time-step 0; 6-12 seconds 属于 time-step 1

		for (int i = 0; i < N_up.size() - 1; i++) {
			for (int j = 0; j < N_down.size() - 1; j++) {

				// TODO changing Param.time units, dt & seconds is huge difference
				double b1 = (wb*uf)* N_up.get(i) + ((wb*uf)*(current_time) - (wb*uf)*(current_time - (N_up.size() - (i+1))))*(N_up.get(i+1) - N_up.get(i));  //  ---> refer 推算 04-10-2022 Page-11

				double b2 = (wb*uf)* N_down.get(j) + ((wb*uf)*(current_time) - uf*L - (wb*uf)*(current_time  - (N_down.size() - (j+1))))*(N_down.get(j+1) - N_down.get(j)) + wb*uf*K_jam*L ; //  ---> refer 推算 04-10-2022 Page-11

				double b = (b1 - b2);

				double a  = (uf*(N_down.get(j+1) - N_down.get(j)) - wb*uf*K_jam + wb*(N_up.get(i+1) - N_up.get(i)));

				// TODO:---> Te puts SolveLinearEquation in Link.java
				temp_xs = SolveLinearEquation(a, b);
				//System.out.println("temp xs is:" + temp_xs); // output 都是NaN ---> 有多少个NaN? 32个 基于 i < N_up2.size() - 1 &  j < N_down2.size() - 1

				double up_stream_time_check;
				double down_stream_time_check;

				up_stream_time_check = current_time - (temp_xs/uf); // ---> (t - xs / uf): track-back time for upstream N(up_stream_time_check);
				down_stream_time_check = current_time - (L - temp_xs) / wb; // ---> (t - (L-x_s)/wb): track-back time for downstream N(down_stream_time_check);

				int up_stream_time_lowerbound;
				int up_stream_time_upperbound;

				// check i， 是否 index 到了正确的值 ---> 看手算笔记 04-10-2022 page-7
				up_stream_time_lowerbound = current_time - (N_up.size() - (i+1));  // i 是 index 这里可能会有bug
				up_stream_time_upperbound = current_time - (N_up.size() - (i+2));

				int down_stream_time_lowerbound;
				int down_stream_time_upperbound;

				// check j， 是否 index 到了正确的值
				down_stream_time_lowerbound = current_time  - (N_down.size() - (j+1)); // j
				down_stream_time_upperbound = current_time  - (N_down.size() - (j+2)); // j

				if (up_stream_time_check >= up_stream_time_lowerbound && up_stream_time_check <= up_stream_time_upperbound && down_stream_time_check >=  down_stream_time_lowerbound && down_stream_time_check <= down_stream_time_upperbound) {

					correct_x_shock.add(temp_xs);

				}

			}

		}

		// x_s2 就是 x_shock
		for (double item: correct_x_shock) {
			//System.out.println("print out each item in correct_x_shock2: " + item);
		}
		HashSet<Double> h = new HashSet<Double>(correct_x_shock);
		LinkedList<Double> remove_duplicated_correct_x_shock = new LinkedList<Double>(h);
		System.out.println("After removing duplicated value in correct_x_shock2, we get single & valid x_shock value at link" +  " " + this.getId() + " " + " which is: " +  " " +remove_duplicated_correct_x_shock + " " + "and current time is" + " " + current_time);

		LinkedList<Double> x_shock = remove_duplicated_correct_x_shock;
		return x_shock;

	}


	//TODO: Te you must changing Param.time logic which is different from MWL code
	public LinkedList<Double> shockN(){

		// TODO changing uf wb units
		double uf = ((this.getFFSpeed()/3600)*(Params.dt)); // we need per/dt here
		//double wb = (this.getBackwardSpeed()/3600)*(Params.dt);
		//double K_jam = Params.JAM_DENSITY;
		// TODO changing length units, Yash is ft & MWL is mile
		//double L = (this.getLength()/5280);

		LinkedList<Double> x_shock = this.ShockWaveDetection();

		// TODO: Param.time should be per/dt
		int current_time = (int) Math.floor(Params.time/Params.dt); // 感觉向下取整比较合理 0-6 seconds 属于 time-step 0; 6-12 seconds 属于 time-step 1
		double trackback_time = (current_time - (x_shock.get(0) / uf)); // 10 - (0.14/0.05) = 7.2;
		System.out.println("trackback_time is at link" + " " + this.getId()+ " " + "is" + " " +  trackback_time + " "  + "and current time is" + " " + current_time);

		LinkedList<Double> x_shock_N = new LinkedList<Double>();
		double temp_x_shock_N = 0;

		int downbound_int_time = (int) Math.floor(trackback_time); // down_int_time = 7
		int upbound_int_time = (int) Math.ceil(trackback_time); // up_int_time = 8

		// ---> 如果 trackback_time 就是 7.00 8.00 呢？？？
		// 找一个手稿上的数值再check
		if(downbound_int_time != upbound_int_time ) {

			int max_index = N_up.size() - 1;
			int min_index = 0;

			// (difference of index) = (difference of time) 差值是 一样的
			// 参考 page-58 中的值
			// TODO: Param.time should be per/dt
			int downbound_int_time_index =  max_index - (current_time - downbound_int_time); // 4 - (10 -7) = 1; index = 1, 对应 time = 7 N_up(t=7) = N_up(index=1) = 14
			int upbound_int_time_index = max_index - (current_time - upbound_int_time); // 4 - （10 - 8）= 4 - 2 = 2; index = 2, 对应  time = 8 N_up(t=8）= N_up(index=

			double dN = N_up.get(upbound_int_time_index) -  N_up.get(downbound_int_time_index);
			double dt = trackback_time - (double)downbound_int_time;

			temp_x_shock_N = (N_up.get(downbound_int_time_index) +  dt * dN); // math.round make 14.4000000 to 14.0000000
			x_shock_N.add(temp_x_shock_N);

		} else {


			int max_index = N_up.size() - 1;
			// if trackback_time 就是一个整数 对应一个
			// TODO: This is place that you make mistake twice
			//  * Te you should return the index not the time (4 - (10 - 6)), otherwise you will "IndexOutOfBoundsException"
			temp_x_shock_N = N_up.get(max_index -(current_time - (int)trackback_time));
			x_shock_N.add(temp_x_shock_N);
		}

		return x_shock_N;

	}



	public LinkedList<Double> downstreamBreakpointDetection() {
		//traffic parameters
		// TODO changing uf wb units, per/dt & per/seconds is huge difference
		// double uf = ((this.getFFSpeed()/3600)*(Params.dt));
		double wb = (this.getBackwardSpeed()/3600)*(Params.dt);
		//double K_jam = Params.JAM_DENSITY;
		// TODO changing length units, Yash is ft & MWL is mile
		double L = (this.getLength()/5280);

		// 变量的初始化
		LinkedList<Double> x_break_downstream = new LinkedList<Double>();
		double temp_backward_position = 0;
		int max_N_down_index = (N_down.size() - 1);

		// ---> N_down flow change point checking
		// temp_backward2 需要与 temp_xs2 比较后 才能看是否能存储下来
		for (int j =0; j < N_down.size() - 2; j++) {
			// TODO: 这个找backward speed breakpoint的 if 条件很重要: 这个 跟 upstream break point 不一样，我们找的是 red light change to green light's backward wave
			if( (N_down.get(j+2) - N_down.get(j+1)) > 0 && (N_down.get(j+1) - N_down.get(j)) == 0 ) {
				//System.out.println("Downstream flow backward speed start point index is:" + (j+1)); // ---> corresponding to LinkedList index not the real time

				temp_backward_position = L - wb*( max_N_down_index - (j+1)); //  L - wb * (dt)
				//System.out.println("Backward flow reach point is:" + temp_backward);

				// LinkedList correct_xs2 add valid value
				x_break_downstream.add(temp_backward_position);
			}

		}

		return  x_break_downstream; // this is a LinkedList

	}

	public LinkedList<Double> backwardBreakN() {

		// traffic parameters;
		// TODO changing uf wb units, per/dt & per/seconds is huge difference
		//double uf = ((this.getFFSpeed()/3600)*(Params.dt));
		double wb = (this.getBackwardSpeed()/3600)*(Params.dt);
		double K_jam = Params.JAM_DENSITY;
		// TODO changing length units, Yash is ft & MWL is mile
		double L = (this.getLength()/5280);

		LinkedList<Double> x_break_downstream = this.downstreamBreakpointDetection();

		LinkedList<Double> x_break_downstream_N = new LinkedList<Double>();
		// trackback to the break point index
		int max_N_down_index = (N_down.size() - 1);
		// x_break_downstream 如果存在 则 只存在一个值 所以取第一个index = 0
		int backwardBreakN_index = (int) (Math.round ((x_break_downstream.get(0)  + wb * max_N_down_index - L)/wb) );  // 正确使用: Math.round（包含所有计算）
		//System.out.println("value in the bracket is: " + ((x_break_downstream.get(0)  + wb * max_N_down_index - L)/wb));
		System.out.println("backwardBreakN_index is: " + backwardBreakN_index + " " + "at linkid " +  " " + this.getId());

		// TODO: Te you make mistakes here twice, be careful of downstream math N(t-(L-x_s)/wb) + K(L-xs)
		//    * you miss this part " K(L-xs) " twice !!!
		double  temp_x_break_downstream_N = N_down.get(backwardBreakN_index) + K_jam * (L - x_break_downstream.get(0) );
		x_break_downstream_N.add(temp_x_break_downstream_N);

		return x_break_downstream_N;

	}


	public LinkedList<Double> upstreamBreakpointDetection(){
		// traffic parameters
		// TODO changing uf wb units, per/dt & per/seconds is huge difference
		double uf = ((this.getFFSpeed()/3600)*(Params.dt));
		//double wb = (this.getBackwardSpeed()/3600)*(Params.dt);
		//double K_jam = Params.JAM_DENSITY;
		// TODO changing length units, Yash is ft & MWL is mile
		//double L = (this.getLength()/5280);

		LinkedList<Double> x_break_upstream = new LinkedList<Double>();
		double temp_freeflow_trackback = 0;
		double dt = 0;

		for (int ii = 0; ii < N_up.size() - 2; ii++) {

			/**
			 * This if condition is different when checking break point for N_down
			 */
			if ((N_up.get(ii+2) - N_up.get(ii+1)) !=  (N_up.get(ii+1) - N_up.get(ii))) {
				//System.out.println("upstream flow change time point index is:" + (ii+1));
				// TODO: dt should be time-step
				dt = ((N_up.size() - 1) - (ii+1));
				temp_freeflow_trackback = uf * dt; // 速度 uf2 * 经历的时间 ((N_up2.size() - 1) - (ii+1))
				x_break_upstream.add(temp_freeflow_trackback);
			}

		}

		return x_break_upstream;
	}

	// 需要一个loop loop出 LinkedList x_break_upstream 里面所有的值 然后反算index
	// TODO: x_break_upstream 里面的解还有顺序问题
	// testing data input ---> refer handwriting draft page-58
	// 这个改成 return double 是否更合适， 如果还是LinkedList 想想怎么code
	public double forwardBreakN(double single_x_break_upstream) {

		//traffic parameters
		// TODO: uf should be per/dt

		double uf = ((this.getFFSpeed()/3600)*(Params.dt));

		// TODO: put x_start in some better classes;
		double x_start = 0.00;

		int max_N_up_index = (N_up.size() - 1);
		double dt = 0;
		double dx = 0;

		double x_break_upstream_N =  0;

		//TODO: 0.15000000000000002 输入 return 值有问题 解决精度问题 先用多余的弄了
		//  * 精度问题
		//  * dt is difference meaning of Params.dt 这里的dt 是用来 track back time的
		dx = single_x_break_upstream - x_start;
		//System.out.println("dx ----> " + dx);
		dt = dx / uf;
		//System.out.println("dt ----> " + dt);

		// track back to index
		// TODO: Possible bug maybe due to precision problems
		int temp_index = (int)(Math.round(max_N_up_index - dt)); // Math.round (5.9999999999999964) = 5.0 Math.round (6.0000000000001) = 6.0
		//System.out.println("(max_N_up_index - dt) is:" + (max_N_up_index - dt));
		System.out.println("temp_index is: " + " " + temp_index + " for linkid is " + " " + this.getId());
		double temp_x_break_upstream_corresponding_N = N_up.get(temp_index);

		x_break_upstream_N = temp_x_break_upstream_corresponding_N;
		return x_break_upstream_N;

	}




	public HashMap<Double, Double> get_X_N_pairs() {

		// traffic params
		//double uf = this.getFFSpeed();
		//double wb = Params.wb;
		//double K_jam = Params.K_jam;
		// TODO: L is mile based on my original get_X_N_pairs()
		double L = (this.getLength()/5280);
		// 从 x_start 到 x_end 顺序 put key & value pairs

		HashMap<Double, Double> X_N_pairs = new HashMap<Double, Double>();

		double x_start = 0;
		double x_end = L; // ---> this should be the length of this link or we say current link, notice this this this this
		int max_N_up_index = N_up.size() - 1;
		//System.out.println(max_N_up_index);
		int max_N_down_index = N_down.size() - 1;;

		double x_start_N = N_up.get(max_N_up_index);
		double x_end_N = N_down.get(max_N_down_index);

		// 顺序1 x_start
		X_N_pairs.put(DoubleE.round(x_start, 6), x_start_N);

		// 顺序2 x_break_upstream
		//
		LinkedList<Double> x_break_upstream = this.upstreamBreakpointDetection();

		if (x_break_upstream.size() >  0) {
			// for each x_break_upstream get corresponding N

			for (int i = 0; i <= x_break_upstream.size() - 1; i++) {
				double single_x_break_upstream = x_break_upstream.get(i);
				System.out.println("we have upstream break boundary point, the value is: " + single_x_break_upstream);

				double temp_x_break_upstream_N = this.forwardBreakN(single_x_break_upstream);
				X_N_pairs.put(DoubleE.round(single_x_break_upstream, 6), temp_x_break_upstream_N);
			}

		} else {
			System.out.println("we do not have upstream break boundary point");
		}


		//顺序3 x_shock
		LinkedList<Double> x_shock = this.ShockWaveDetection();
		System.out.println("we have x_shock solution, value is" + x_shock);
		// check 是否有 shockwave
		// TODO: x_shock solution 可能会是 link length 的长度 AKA 的 link 的 endpoint 这个点要 remove ？？？？ 但是有没有更好的的方法 下次开会讨论
		//    * 算 simulation reach time-step = 15 link b case的时候 x_break_upstream 可能会 > x_shock 这是不符合实际的 也应该剔除
		//    * LinkedList class size() method returns the number of elements present in the given list. If a given LinkedList is empty then its size will be 0(zero).
		if (x_shock.size() > 0 && x_shock.get(0)!= x_end ) {
			LinkedList<Double> x_shock_N = shockN();
			System.out.println("We have x_shock solution, which means we face congestion, and x_shock position is: " + x_shock + " and correspoding # of vehicles are " + x_shock_N + "!!!!! and we will test whether we have backward boundary" + "now we are at link" + this.getId());
			// x_shock is a linkedlist we need loop them to get each value
			// However, since at each time for each link we have only one x_shock and x_shock_N



			X_N_pairs.put( DoubleE.round(x_shock.get(0), 6 ), x_shock_N.get(0));

			// 想确认是否有 x_break_upstream 的值
			if (x_break_upstream.size() > 0) {
				for (int i = 0; i <= x_break_upstream.size() - 1; i++) {
					double single_x_break_upstream = x_break_upstream.get(i);
					System.out.println("we have upstream break boundary point, since we also have x_shock point, we need to compare them. If x_break_upstream > x_shock, we need to remove x_break_upstream " +  "now we are at link" + this.getId());
					if (single_x_break_upstream >  x_shock.get(0)) {
						System.out.println(" Yes !!!!! We find some in-valid x_break_upstream point" + " , now we are at link" + this.getId());
						// have map remove key-value pairs
						X_N_pairs.remove(single_x_break_upstream);
					}

				}

			}


		} else {
			System.out.println("we do not have shockwave, so we don't have backward boundary either && We do not need to check x_shock and and x_break_upstream");
		}



		// 顺序4 x_break_downstream
		// check 是否有 shockwave 和 backward boundary 同时存在
		// if (backward boundary > shockwave) save both
		// Otherwise, only save shockwave
		LinkedList<Double> x_break_downstream = this.downstreamBreakpointDetection();
		System.out.println("we have x_break_downstream solution, value is" + " " + x_break_downstream + " " +  "now we are at link" + this.getId());
		// Notice --->  x_break_downstream value > x_shock value 才认为 这个 x_break_downstream 是有效值
		if (x_shock.size() > 0 && x_break_downstream.size() > 0 && x_break_downstream.get(0) > x_shock.get(0)) {
			LinkedList<Double> x_break_downstream_N =  this.backwardBreakN() ;
			System.out.println("We have valid x_shock solution and x_break_downstream solution, so we will save them" + ", now we are at link" + " " + this.getId());
			X_N_pairs.put(DoubleE.round(x_break_downstream.get(0), 6), x_break_downstream_N.get(0));
		} else {
			if (x_break_downstream.size() > 0) {
				System.out.println("This is an invalid x_break_downstream solution since x_break_donwstream < x_shock " + ", now we are at link" + " " + this.getId());
			} else {
				System.out.println("We don't have backward boundary AKA !downstreamBreakpoint! now we are at link" + " " +  this.getId());
			}
		}


		//顺序5: x_end
		X_N_pairs.put(DoubleE.round(x_end,6), x_end_N);

		return sortbykey(X_N_pairs);

	}


	public static HashMap<Double, Double> sortbykey(HashMap<Double, Double> map)
	{
		ArrayList<Double> sortedKeys = new ArrayList<Double>(map.keySet());
		HashMap<Double, Double> sorted_X_N_pairs = new HashMap<>();

		// 上面是 x_end: 0.2 mile 下面是 x_start: 0.0 mile
		// reverse 是反转 不是降序
		// TODO: check this 易错
		Collections.sort(sortedKeys, Collections.reverseOrder());
		// Collections.sort(sortedKeys);

		for (double x : sortedKeys) {
			sorted_X_N_pairs.put(x, map.get(x));
			System.out.println("Key ======= " + x + ", Value ======= " + map.get(x));
		}

		return sorted_X_N_pairs;

	}


	/** This is My APWBP code
	 */
	public static double calculatePressure_upSource(EntryLink upLink, CTM downLink, double turningProportion) {

		//double EntryLink_K = (Params.K_c)/2;

		// TODO: L should be mile
		double downLink_length = (downLink.getLength()/5280);
		HashMap<Double, Double> sorted_X_N_pairs_downLink = downLink.get_X_N_pairs();

		// for upstream is EntryLink, which is the Source Link in Jabari's paper
		// TODO: Occupancy() is number of vehicle on links; This is Jabari's paper pressure calculation when upLink is Src or EntryLink
		//  * checking the math for this pressure calculation ---> when a is Asrc a to b has turningProportion ????


		// TODO ----------------------------------------- Checking this line -----------------------------------------
		//    * 根据Varaiya 的公式来看 可能根本就没有这个值 可能就是Occupancy ！！! 到底有没有 turningProportion 呢？？？
		// double upLink_part_value = upLink.getOccupancy()  * Params.Constant_C; // upSrc to downstream turningProportion
		double upLink_part_value = upLink.getOccupancy()  * Params.Constant_C * turningProportion;
		// System.out.println("upLink getOccupancy is ==>" + " " + upLink.getOccupancy() + " "  + " Constant C is " + " " + Params.Constant_C + "" + "upLink a (EntryLink or Source Node)");

		// for downstream
		// for downLink
		ArrayList<Double> sortedKeys_downLink = new ArrayList<Double>(sorted_X_N_pairs_downLink.keySet());
		Collections.sort(sortedKeys_downLink);
		//Collections.sort(sortedKeys, Collections.reverseOrder());
		List<Double> indexes_downLink = new ArrayList<Double>(sortedKeys_downLink); // --> Notice: List of sorted_X_N_pairs indexes is: [0.0, 0.15, 0.1, 0.2] if we use Collections.sort(sortedKeys);
		System.out.println("sorted X_N_pairs indexes for downstream Link (link_b) in jabari's paper, when link_a is entry link) List is >>>>>> " + " " + indexes_downLink + " " + " and the link id is" + " " + downLink.getId());

		double temp_Density_downLink;
		double dN_downLink;
		double dX_downLink;
		double temp_Intergration_downLink1 = 0;
		double temp_Intergration_downLink2 = 0;
		double downLink_part_value = 0;


		for (int i = 0; i < sorted_X_N_pairs_downLink.size() - 1; i++) {
			dN_downLink = sorted_X_N_pairs_downLink.get(indexes_downLink.get(i)) -  sorted_X_N_pairs_downLink.get(indexes_downLink.get(i+1));
			dX_downLink = indexes_downLink.get(i+1) - indexes_downLink.get(i);
			temp_Density_downLink = dN_downLink/dX_downLink;
			//System.out.println("temp_Density_downLink from x =  " + " " +  indexes_downLink.get(i) + " " + "-----to----" +  indexes_downLink.get(i+1) + " " + " ======> " + temp_Density_downLink );
			temp_Intergration_downLink1 += temp_Density_downLink * (indexes_downLink.get(i+1) - indexes_downLink.get(i));
			//System.out.println("temp_Intergration_downLink1 ======> " + temp_Intergration_downLink1 );
			//System.out.println("temp_Density_downLink * (1 / downLink_length)======> " +temp_Density_downLink * (1 / downLink_length) );
			//System.out.println("((0.5)* Math.pow(indexes_downLink.get(i+1), 2) - (0.5)* Math.pow(indexes_downLink.get(i), 2)) ======> " + ((0.5)* Math.pow(indexes_downLink.get(i+1), 2) - (0.5)* Math.pow(indexes_downLink.get(i), 2)));
			temp_Intergration_downLink2 += temp_Density_downLink * (1 / downLink_length) * ((0.5)* Math.pow(indexes_downLink.get(i+1), 2) - (0.5)* Math.pow(indexes_downLink.get(i), 2));
			//System.out.println("temp_Intergration_downLink2 ======> " + temp_Intergration_downLink2);

		}


		double temp_Intergration_downLink = temp_Intergration_downLink1 - temp_Intergration_downLink2;

		VehIntersection neigh = (VehIntersection) downLink.getDest();
		Set<TurningMovement> nextTurns =  neigh.getVehicleTurns();
		//System.out.println("\t\t downstream turns: " + nextTurns);
		for (TurningMovement nextTurn: nextTurns){
			if(nextTurn.getIncomingLink().equals(downLink)){
				TurningMovement downstream_turn = nextTurn;
				//System.out.println("turningProportion for downLink b to its down stream link c is =====>" + downstream_turn.getTurningProportion() );
				// ---> 积分后的值乘以 b 到所有 c的 turningProportion
				//downLink_part_value += downstream_turn.getTurningProportion() * downLink_part_value;
				downLink_part_value +=  Params.Constant_C * (downstream_turn.getTurningProportion() * downstream_turn.getTurningProportion()) * temp_Intergration_downLink;
				//System.out.println("I am checking the every for loop for downstream part value for pressure calculation in Jabari's paper is =====>" + downLink_part_value );

			}
		}

		double pressure_Value = Math.abs(upLink_part_value - downLink_part_value);
		System.out.println("The pressure value for upLink a calculated from my paper APWBP is ---------->" + upLink_part_value);
		System.out.println("The pressure value for downLink a calculated from my paper APWBP is ---------->" + downLink_part_value);
		System.out.println("The pressure value calculated from my paper APWBP is (when upstream is entrylink) ---------->" + " " + pressure_Value + " " + "and the uplink id is" + " " + upLink.getId() + "and the downlink id is" + " " + downLink.getId());
		return pressure_Value;

	}





	// public static double calculatePressure_normal(LTM upLink, LTM downLink, double turningProportion, int cur_time, LinkedList<Double> try_N_up, LinkedList<Double> try_N_down)
	public double calculatePressure_normal(CTM upLink, CTM downLink, double turningProportion) {


		// ---> Be care for getLength() since upLink length is difference from downLink
		// TODO: be careful L's units
		double upLink_length = (upLink.getLength()/5280); // change this to upLink.getLength() when simulation start
		double downLink_length = (downLink.getLength()/5280); // change this to downLink.getLength() when simulation start

		//int linkDensitysize = sorted_X_N_pairs.size() - 1;
		//ArrayList<Double> linkDensity = new ArrayList<>(linkDensitysize);

		HashMap<Double, Double> sorted_X_N_pairs_upLink = upLink.get_X_N_pairs();
		HashMap<Double, Double> sorted_X_N_pairs_downLink = downLink.get_X_N_pairs();

		/** For upLink
		 */
		ArrayList<Double> sortedKeys_upLink = new ArrayList<Double>(sorted_X_N_pairs_upLink.keySet());
		Collections.sort(sortedKeys_upLink);
		//Collections.sort(sortedKeys, Collections.reverseOrder());
		List<Double> indexes_upLink = new ArrayList<Double>(sortedKeys_upLink); // --> Notice: List of sorted_X_N_pairs indexes is: [0.0, 0.15, 0.1, 0.2] if we use Collections.sort(sortedKeys);
		System.out.println("sorted X_N_pairs indexes for upstream Link (link_a in jabari's paper) List is ======>  " + " " + indexes_upLink + "and the uplink id is" + " " + upLink.getId());

		double temp_Density_upLink;
		double dN_upLink;
		double dX_upLink;
		double temp_Intergration_upLink;
		double upLink_part_value = 0;
		// Math.pow(a,b) -> a的b次方 不是 Math.sqrt() !!!!!
		for (int i = 0; i < sorted_X_N_pairs_upLink.size() - 1; i++) {
			dN_upLink = sorted_X_N_pairs_upLink.get(indexes_upLink.get(i)) -  sorted_X_N_pairs_upLink.get(indexes_upLink.get(i+1)); // i < 4 ; i = 3 i+1 = 4 这个是index 是对的
			dX_upLink = indexes_upLink.get(i+1) - indexes_upLink.get(i);
			temp_Density_upLink = dN_upLink/dX_upLink;
			//System.out.println("temp_Density_upLink from x =  " + " " +  indexes_upLink.get(i) + " " + "-----to----" +  indexes_upLink.get(i+1) + " " + " ======> " + temp_Density_upLink );
			//System.out.println(" x ======> " + indexes_upLink.get(i+1));
			//System.out.println(" x ======> " + indexes_upLink.get(i));
			// System.out.println(" x^2 ======> " + Math.pow(indexes_upLink.get(i+1), 2) );
			// System.out.println(" x^2 ======> " + Math.pow(indexes_upLink.get(i), 2) );

			temp_Intergration_upLink = (0.5)* Math.pow(indexes_upLink.get(i+1), 2) - (0.5)* Math.pow(indexes_upLink.get(i), 2);
			//  System.out.println("temp_Intergration_upLink ======> " + temp_Intergration_upLink);
			//linkDensity.add(temp_Density_upLink);
			//  System.out.println("temp_Density_upLink * temp_Intergration_upLink ======> " + (temp_Density_upLink * temp_Intergration_upLink));
			// TODO: Constant_C may be a trick in paper, I am not sure how it will influence results, Maybe need to ask Michael. ---> Michael said just setting to 1
			//  System.out.println("C_ab * (1/la) ======> " + (Params.Constant_C)*(1/upLink_length));
			upLink_part_value += (Params.Constant_C)*(1/upLink_length)*(temp_Density_upLink * temp_Intergration_upLink);

			//System.out.println("upLink_part_value ======> " + " " + upLink_part_value + " " + " link id is" + " " + upLink.getId() );
		}


		upLink_part_value = upLink_part_value * turningProportion;


		/** For downLink
		 */
		ArrayList<Double> sortedKeys_downLink = new ArrayList<Double>(sorted_X_N_pairs_downLink.keySet());
		Collections.sort(sortedKeys_downLink);
		//Collections.sort(sortedKeys, Collections.reverseOrder());
		List<Double> indexes_downLink = new ArrayList<Double>(sortedKeys_downLink); // --> Notice: List of sorted_X_N_pairs indexes is: [0.0, 0.15, 0.1, 0.2] if we use Collections.sort(sortedKeys);
		System.out.println("sorted X_N_pairs indexes for downstream Link (link_b in jabari's paper) List is ======> " + " " + indexes_downLink + "and the downlink id is" + " " + downLink.getId());


		double temp_Density_downLink;
		double dN_downLink;
		double dX_downLink;
		double temp_Intergration_downLink1 = 0;
		double temp_Intergration_downLink2 = 0;
		double downLink_part_value = 0;


		for (int i = 0; i < sorted_X_N_pairs_downLink.size() - 1; i++) {
			dN_downLink = sorted_X_N_pairs_downLink.get(indexes_downLink.get(i)) -  sorted_X_N_pairs_downLink.get(indexes_downLink.get(i+1)); ;
			dX_downLink = indexes_downLink.get(i+1) - indexes_downLink.get(i);
			// TODO: please make sure temp_Density units, should be vehs per mile
			temp_Density_downLink = dN_downLink/dX_downLink;
			//System.out.println("temp_Density_downLink from x =  " + " " +  indexes_downLink.get(i) + " " + "-----to----" +  indexes_downLink.get(i+1) + " " + " ======> " + temp_Density_downLink );
			// TODO: This is just occupancy, number of vehs on the link
			temp_Intergration_downLink1 += temp_Density_downLink * (indexes_downLink.get(i+1) - indexes_downLink.get(i));

			assert Double.compare(temp_Intergration_downLink1, getOccupancy()) == 0;

			//System.out.println("temp_Intergration_downLink1 ======> " + temp_Intergration_downLink1 );
			//System.out.println("temp_Density_downLink * (1 / downLink_length)======> " +temp_Density_downLink * (1 / downLink_length) );
			//System.out.println("((0.5)* Math.pow(indexes_downLink.get(i+1), 2) - (0.5)* Math.pow(indexes_downLink.get(i), 2)) ======> " + ((0.5)* Math.pow(indexes_downLink.get(i+1), 2) - (0.5)* Math.pow(indexes_downLink.get(i), 2)));
			temp_Intergration_downLink2 +=  temp_Density_downLink * (1 / downLink_length) * ((0.5)* Math.pow(indexes_downLink.get(i+1), 2) - (0.5)* Math.pow(indexes_downLink.get(i), 2));
			//System.out.println("temp_Intergration_downLink2 ======> " + temp_Intergration_downLink2);


		}

		double temp_Intergration_downLink = temp_Intergration_downLink1 - temp_Intergration_downLink2;

		// TODO:                                                                                                     ^
		// TODO:                                                                                                    /
		// TODO:                                                                                                  c2
		// TODO:                                                                                                 /
		// TODO: checking the math for pressure calculation: ----a----> you are controlling this node ----b----> ----c1---->

		VehIntersection neigh = (VehIntersection) downLink.getDest();
		Set<TurningMovement> nextTurns =  neigh.getVehicleTurns();
		//System.out.println("\t\t downstream turns: " + nextTurns);
		for (TurningMovement nextTurn: nextTurns){
			if(nextTurn.getIncomingLink().equals(downLink)){
				TurningMovement downstream_turn = nextTurn;
				//System.out.println("turningProportion for downLink b to its down stream link c is =====>" + downstream_turn.getTurningProportion() );
				// ---> 积分后的值乘以 b 到所有 c的 turningProportion
				//downLink_part_value += downstream_turn.getTurningProportion() * downLink_part_value;
				downLink_part_value +=  Params.Constant_C * (downstream_turn.getTurningProportion() * downstream_turn.getTurningProportion()) * temp_Intergration_downLink;
				//System.out.println("I am checking the every for loop for downstream part value for pressure calculation in Jabari's paper is =====>" + downLink_part_value );

			}
		}

		System.out.println("The pressure value for upLink a calculated from my paper APWBP is ---------->" + upLink_part_value);
		System.out.println("The pressure value for downLink a calculated from my paper APWBP is ---------->" + downLink_part_value);

		double pressure_Value = Math.abs(upLink_part_value - downLink_part_value);
		System.out.println("The pressure value calculated from my paper APWBP is (when upstream is normal link) ---------->" + " " + pressure_Value + " " + "and the uplink id is" + " " + upLink.getId() + " " + " and the downlink id is" + " " + downLink.getId());

		return pressure_Value ;

	}












	// --------------------------------------------------- Jabari's Algorithm ----------------------------------------//

	// TODO "x_start" --- cell index= 0 ---> "x_cell0" --- cell index = 1 ---> "x_cell1" --- cell index = 2 ---> "x_cell2" --- cell index = 3 ---> "x_cell3" --- cell index = 4 ---> "x_cell4"
	//    * We use n not N, because n is the number of vehs in Cell, not cumulative number of vehs in cell
	public List<Double> get_Cell_X_List(){
		int number_of_cells = cells.length; // return number of cells in current CTM link

		double x_start = 0;

		// TODO: be care of getCellLength()'s units
		//    *: 这里很容易出错！！！
		double dx = this.getCellLength(); // return length (mi) of each cell in CTM-Link
		List<Double> Cell_X_List = new ArrayList<Double>();

		/**
		 * threeLinks example with L = 0.2mile
		 * x_start = 0;
		 * cell 0 = [x_start   # of vehs   x_cell_0]  x_cell_0 = x_start + dx
		 * cell 1 = [x_cell_0  # of vehs   x_cell_1]  x_cell_1 = x_start + dx
		 * cell 2 = [x_cell_1  # 0f vehs   x_cell_2]  x_cell_2 = x_start + dx
		 * cell 3 = [x_cell_2  # of vehs   x_cell_3]  x_cell_3 = x_start + dx
		 *
		 * Cell_X_List = [0, 0 + dx, 0 + 2dx, 0 + 3dx, 0 + 4dx] size = 5
		 */

		// add a size checking logic
		Cell_X_List.add(x_start); // index = 0, put in x_start
		double temp_cell_x = 0;
		// start form index = 1
		/** (index = 0,x_start) --> (index = 1, x_start + dx)"first cell end"  --> (index = 2, x_start + 2dx) "second cell end" --> (index = 3, x_start + 3dx) "third cell end" --> (index = 4, x_start + 4dx) "forth cell end"
		 * i <= number_of_cells ---> i <= 4
		 */
		for (int i = 1; i <= number_of_cells; i ++){
			temp_cell_x += dx;
			Cell_X_List.add(temp_cell_x); // index = 1, put in 0 + dx; index = 2, put in 0 + dx + dx;
		}

		// cells.length is number of cells for this CTM-Link
		// cells_x
		// Cell_X_List size should be 1 more than cell.length
		if (Cell_X_List.size() !=  cells.length + 1){
			System.out.println("Cell_X_List size is not correct");
		}
		System.out.println("Cell_X_List for current link id"  + " " + this.getId() + " " + Cell_X_List);
		return Cell_X_List;
	}


	/**
	 * @return current CTM-LINK cells
	 */
	public Cell[] getCells() {
		return cells;
	}


	// TODO 面向对象编程 区别upLink downLink cells Cell 别index 错了！！！
	public static double calculatePressure_upSource_Jabari(EntryLink upLink, CTM downLink, double turningProportion) {

		// TODO: L should be mile
		double downLink_length = (downLink.getLength()/5280);
		// checking
		System.out.println("downLink_length check ---> " + " " + downLink_length);
		// TODO: Be care for the units of below
		double downLink_eachCell_length  = downLink.getCellLength();
		// checking
		System.out.println("downLink_eachCell_length check ---> " + " " + downLink_eachCell_length);

		// for upstream
		double upLink_part_value = upLink.getOccupancy()  * Params.Constant_C * turningProportion;
		// checking
		System.out.println("upLink getOccupancy is ==>" + " " + upLink.getOccupancy() + " "  + " Constant C is " + " " + Params.Constant_C + "" + "upLink a (EntryLink or Source Node)");


		// for downstream
		double temp_Density_downLink; // this is for each cell in CTM ===> # of vehs / dx (getCellLength())
		double temp_Intergration_downLink1 = 0;
		double temp_Intergration_downLink2 = 0;
		double downLink_part_value = 0;
		List<Double> Cell_X_index_downLink = downLink.get_Cell_X_List(); // return [0, 0 + dx, 0 + 2dx, 0 + 3dx, 0 + 4dx]
		// checking
		System.out.println("downLink's number of cells in total ---> " + " " + downLink.getCells().length);

		for (int cell_index = 0; cell_index < downLink.getCells().length; cell_index ++){
			temp_Density_downLink = downLink.getCells()[cell_index].getOccupancy() / downLink_eachCell_length;

			temp_Intergration_downLink1 += temp_Density_downLink * (downLink.getCellLength() *(cell_index + 1) - downLink.getCellLength() *(cell_index));
			temp_Intergration_downLink2 +=  temp_Density_downLink * (1 / downLink_length) * ((0.5)* Math.pow(downLink.getCellLength() * (cell_index + 1), 2) - (0.5)* Math.pow(downLink.getCellLength() *(cell_index), 2));
		}


		double temp_Intergration_downLink = temp_Intergration_downLink1 - temp_Intergration_downLink2;


		VehIntersection neigh = (VehIntersection) downLink.getDest();
		Set<TurningMovement> nextTurns =  neigh.getVehicleTurns();
		//System.out.println("\t\t downstream turns: " + nextTurns);
		for (TurningMovement nextTurn: nextTurns){
			if(nextTurn.getIncomingLink().equals(downLink)){
				TurningMovement downstream_turn = nextTurn;
				//System.out.println("turningProportion for downLink b to its down stream link c is =====>" + downstream_turn.getTurningProportion() );
				// ---> 积分后的值乘以 b 到所有 c的 turningProportion
				//downLink_part_value += downstream_turn.getTurningProportion() * downLink_part_value;
				downLink_part_value +=  Params.Constant_C * (downstream_turn.getTurningProportion() * downstream_turn.getTurningProportion()) * temp_Intergration_downLink;
				//System.out.println("I am checking the every for loop for downstream part value for pressure calculation based on Jabari's Algo is =====>" + downLink_part_value );
			}
		}


		System.out.println("The pressure value for upLink a calculated based on Jabari's Algo is ---------->" + " " +  upLink_part_value);
		System.out.println("The pressure value for downLink a calculated based on Jabari's Algo is ---------->" + " " +  downLink_part_value);

		double pressure_Value = Math.abs(upLink_part_value - downLink_part_value);
		System.out.println("The pressure value calculated based on Jabari's Algo is (when upstream is normal link) ---------->" + " " + pressure_Value + " " + "and the uplink id is" + " " + upLink.getId() + " " + " and the downlink id is" + " " + downLink.getId());
		return pressure_Value;


	}


	// TODO 面向对象编程 区别upLink downLink cells Cell 别index 错了！！！
	public double calculatePressure_normal_Jabari(CTM upLink, CTM downLink, double turningProportion) {

		// TODO: be careful L's units
		double upLink_length = (upLink.getLength()/5280); // change this to upLink.getLength() when simulation start
		double downLink_length = (downLink.getLength()/5280); // change this to downLink.getLength() when simulation start

		// TODO: Be care for the units of below
		//    * Because if upstream and downstream ffspd is differenct, each Cell's length is different
		double upLink_eachCell_length = upLink.getCellLength();
		double downLink_eachCell_length  = downLink.getCellLength();

		/**
		 * This for upLink intergration through each cell
		 */
		double temp_Density_upLink; // this is for each cell in CTM ===> # of vehs / dx (getCellLength())
		double temp_Intergration_upLink;
		double upLink_part_value = 0;
		List<Double> Cell_X_index_upLink = upLink.get_Cell_X_List(); // return [0, 0 + dx, 0 + 2dx, 0 + 3dx, 0 + 4dx]

		// cells.length is number of cells for this CTM-Link; eg = 4, last cell index is 4 - 1 =3
		// cell_index < 4
		for (int cell_index = 0; cell_index < upLink.getCells().length; cell_index ++){
			// TODO: please make sure temp_Density units, should be vehs per mile
			temp_Density_upLink = upLink.getCells()[cell_index].getOccupancy() / upLink_eachCell_length;
			System.out.println("cell occupancy in Jabari's Algo from x =  " + " " +  Cell_X_index_upLink.get(cell_index) + " " + "-----to----" +  Cell_X_index_upLink.get(cell_index+1) + " " + " ======> " + upLink.getCells()[cell_index].getOccupancy());
			System.out.println("temp_Density_upLink in Jabari's Algo from x =  " + " " +  Cell_X_index_upLink.get(cell_index) + " " + "-----to----" +  Cell_X_index_upLink.get(cell_index+1) + " " + " ======> " + temp_Density_upLink );
			//
			temp_Intergration_upLink =(0.5)*Math.pow(upLink.getCellLength() * (cell_index + 1), 2) - (0.5)*Math.pow(upLink.getCellLength()*(cell_index), 2);

			// for each cell we have density * temp_Intergration_upLink ---> 这个就是积分
			upLink_part_value +=  (Params.Constant_C)*(1/upLink_length)*(temp_Density_upLink * temp_Intergration_upLink);
		}

		upLink_part_value = upLink_part_value * turningProportion;

		/**
		 * This for downLink intergration through each cell
		 */
		// TODO: please make sure temp_Density units, should be vehs per mile
		double temp_Density_downLink; // this is for each cell in CTM ===> # of vehs / dx (getCellLength())
		double temp_Intergration_downLink1 = 0;
		double temp_Intergration_downLink2 = 0;
		double downLink_part_value = 0;
		List<Double> Cell_X_index_downLink = downLink.get_Cell_X_List(); // return [0, 0 + dx, 0 + 2dx, 0 + 3dx, 0 + 4dx]

		/** 注意 cell */
		for (int cell_index = 0; cell_index < downLink.getCells().length; cell_index ++){
			temp_Density_downLink = downLink.getCells()[cell_index].getOccupancy() / downLink_eachCell_length;
			temp_Intergration_downLink1 += temp_Density_downLink * (downLink.getCellLength() *(cell_index + 1) - downLink.getCellLength() * (cell_index));
			temp_Intergration_downLink2 +=  temp_Density_downLink * (1 / downLink_length) * ((0.5)* Math.pow(downLink.getCellLength() *(cell_index + 1), 2) - (0.5)* Math.pow(downLink.getCellLength() * (cell_index), 2));
		}


		double temp_Intergration_downLink = temp_Intergration_downLink1 - temp_Intergration_downLink2;

		VehIntersection neigh = (VehIntersection) downLink.getDest();
		Set<TurningMovement> nextTurns =  neigh.getVehicleTurns();
		//System.out.println("\t\t downstream turns: " + nextTurns);
		for (TurningMovement nextTurn: nextTurns){
			if(nextTurn.getIncomingLink().equals(downLink)){
				TurningMovement downstream_turn = nextTurn;
				//System.out.println("turningProportion for downLink b to its down stream link c is =====>" + downstream_turn.getTurningProportion() );
				// ---> 积分后的值乘以 b 到所有 c的 turningProportion
				//downLink_part_value += downstream_turn.getTurningProportion() * downLink_part_value;
				downLink_part_value +=  Params.Constant_C * (downstream_turn.getTurningProportion() * downstream_turn.getTurningProportion()) * temp_Intergration_downLink;
				//System.out.println("I am checking the every for loop for downstream part value for pressure calculation based on Jabari's Algo is =====>" + downLink_part_value );
			}
		}


		System.out.println("The pressure value for upLink a calculated based on Jabari's Algo is ---------->" + upLink_part_value);
		System.out.println("The pressure value for downLink a calculated based on Jabari's Algo is ---------->" + downLink_part_value);

		double pressure_Value = Math.abs(upLink_part_value - downLink_part_value);
		System.out.println("The pressure value calculated based on Jabari's Algo is (when upstream is normal link) ---------->" + " " + pressure_Value + " " + "and the uplink id is" + " " + upLink.getId() + " " + " and the downlink id is" + " " + downLink.getId());
		return pressure_Value;
	}


	/**
	 * In Yash new version code, getQueueLength() = getN() * turningProps
	 * In my code, queuelength is total number of vehs waiting to leave at next time
	 * I use turningProps later
	 */
	public double getQueueLength() {
		return getN();
	}

	// number of vehs on the link
	/**
	 * In LTM & CTM paper getN() is the # of vehs that can be send for next time-step
	 * so, it should be sending flow
	 * Yash's sending is different
	 */
	public double getN(){
		return  cells[cells.length-1].getSendingFlow();
	}


	/**
	 * NOTE!!! We should use CTM Link now !!!
	 */
	@Override
	public double getPressure(Link downstreamLink, double turningProportion) {
		// TODO: add case when this link is EntryLink and downstreamLink is LTM
		if (downstreamLink instanceof CTM) {
			CTM dl = (CTM) downstreamLink;
			// TODO: Set params C_ab, C_bc, Pi_bc
			return calculatePressure_normal(this, dl, turningProportion);

		}

		// TODO: what is the other case??
		if (downstreamLink instanceof ExitLink) {
			return 0.0;
		}
		else {
			// TODO: throw error
			new IllegalArgumentException("need downstream to be an LTM link or Exit Link");
		}
		return 0.0;
	}



	/**
	 * NOTE!!! We should use CTM Link now !!!
	 */
	public double getPressure_Jabari(Link downstreamLink, double turningProportion) {
		// TODO: add case when this link is EntryLink and downstreamLink is CTM
		if (downstreamLink instanceof CTM) {
			CTM dl = (CTM) downstreamLink;
			// TODO: Set params C_ab, C_bc, Pi_bc

			return calculatePressure_normal_Jabari(this, dl, turningProportion);

		}

		// TODO: what is the other case??
		if (downstreamLink instanceof ExitLink) {
			return 0.0;
		}
		else {
			// TODO: throw error
			new IllegalArgumentException("need downstream to be an LTM link or Exit Link");
		}

		return 0.0;
	}



}
