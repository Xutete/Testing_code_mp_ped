/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ped;

import java.util.LinkedList;

/**
 * This class propagates flow according to the point queue model (no spatial constraints).
 *
 * @author: Yashveer Bika, Michael W. Levin
 */
public class EntryLink extends Link
{

    private double occupancy;

    public EntryLink() {
        super(0, "entry", null, null, 0, 0, Integer.MAX_VALUE, 1);
    }

    public EntryLink(int id, Node dummySrc, Node dest)
    {
        super(id, "entry", dummySrc, dest, 0, 0, Integer.MAX_VALUE, 1);
        dest.setEntryLink(this);
//        queue = new LinkedList<Vehicle>();

    }


    @Override
    public double getQueueLength() {
        return occupancy;
    }


    //    public EntryLink(int id, Node n)
//    {
//        super(id, null, n, 0, 0, 0, 1);
//
////        queue = new LinkedList<Vehicle>();
//
//    }


    public void reset()
    {
        // TODO:
        // fill this in
//        queue.clear();
    }

    public double getOccupancy()
    {
        // fill this in
        // TODO:
//        return queue.size();
        return occupancy;
    }


    public double getN(){
        return occupancy;
    }


    public void step()
    {
        // fill this in
    }

    public void update()
    {
        // fill this in
    }

    public double getSendingFlow()
    {
        // fill this in
        return occupancy;
    }

    public double getReceivingFlow()
    {
        // fill this in
        return Double.MAX_VALUE;
    }

    public void addFlow(double y)
    {
        occupancy += y;
    }



    @Override
    // TurningProps 对于upSrc Link 还是用的上的
    /*
    public double getPressure(Link downstreamLink, double turningProportion) {
        if (downstreamLink instanceof CTM) {
            CTM dl = (CTM) downstreamLink;
            return CTM.calculatePressure_upSource(this, dl, turningProportion);
            //return LTM.calculatePressure_upSource(this, dl, turningProportion);
        }

        return 0; // What this mean ? never make EntryLink ------> ExitLink directly ????
    }

     */

    // TODO: change to this when runing LTM
    public double getPressure(Link downstreamLink, double turningProportion) {
        if (downstreamLink instanceof LTM) {
            LTM dl = (LTM) downstreamLink;
            return LTM.calculatePressure_upSource(this, dl, turningProportion);
            //return LTM.calculatePressure_upSource(this, dl, turningProportion);
        }

        return 0; // What this mean ? never make EntryLink ------> ExitLink directly ????
    }


    public double getPressure_Jabari(Link downstreamLink, double turningProportion) {
        if (downstreamLink instanceof CTM) {
            CTM dl = (CTM) downstreamLink;
            return CTM.calculatePressure_upSource_Jabari(this, dl, turningProportion);
            //return LTM.calculatePressure_upSource(this, dl, turningProportion);
        }

        return 0; // What this mean ? never make EntryLink ------> ExitLink directly ????
    }



    //TODO: Te add function;
   // public double calculatePressure_normal(LTM upstreamLink, LTM downstreamLink, double turningProportion) {return 0;}

    public void removeFlow(double y)
    {
        occupancy -= y;
    }

    public boolean isEntry() {
        return true;
    }


}