package examples.jmarkov;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.NotUnichainException;

/** Builds two Queues (one VIP and another one for OTHER users)
 * 
 * lambda[]
 *            Arrival rates for VIP and OTHER users to the bank
 * mu[]
 *            Service rate of tellers for each transaction type
 * maxCapacity
 *            Queues capacity
 * M
 *            Number of tellers available in the bank for VIP clients
 * S
 *            Number of tellers available in the bank for OTHER users
 * K   
 *            Different transaction types
 */
        class BankState extends PropertiesState {

                int M; // number of servers serving VIP clients (VIP servers can serve normal clients if there aren't VIP clients waiting to be served )
                int S; // number of servers serving OTHER users (Normal servers are not allowed to serve VIP clients)(Total number of servers (N) minus those serving VIP clients (M))
                int K;// different transaction types                
                
                /** 
                 * Constructs a state for an empty system with M and S servers given K types of transactions.
                 * parameters: M, S, K;
                 */
                
                BankState(int M, int S, int K) {
                        this(new int[K+K+2], 0, 0,K,M,S);
                }
                /**
                 * We identify each State with a vector that counts the 
                 * number of transactions taking place in the VIP servers, then in the OTHER servers, and finally the
                 * number in VIPQueue (in position 2K) and the number in OTHERQueue (in position 2K+1)
                **/
                
                BankState(int[] status, int VIPQsize, int OTHERQsize, int K, int M, int S) {
                        super(K+K+2);
                        this.M = M;
                        this.S = S;
                        this.K = K;
                                
                        for (int j = 0; j < K+K+2; j++) {
                                prop[j] = status[j];}
                        
                        
                        prop[K+K] = VIPQsize;
                        prop[K+K+1] = OTHERQsize;
                        
                }
             

                @Override
                public void computeMOPs(MarkovProcess mp) {
                                        
                        setMOP(mp,"VIPQueue Length", getVIPQSize());
                        setMOP(mp,"OTHERQueue Length", getOTHERQSize());
                
                }
             

                /**
                 * @param i is position in the vector
                 * @return an integer that shows the current state
                 *
                 * 
                 */
                public int getProp(int i) {
                        return prop[i];
                }
               
                /**
                 * Returns the size of the VIP queue 
                 * @return Status of the size of the queue
                 */
                
                public int getVIPQSize() {
                        return prop[K+K];
                }
                  /**
                   * Returns the size of the OTHER queue 
                   * @return Status of the size of the queue
                   */
                
                public int getOTHERQSize() {
                        return prop[K+K+1];
                }
                
                /**
                 * @return  all VIP transactions taking place
                 */
                
                
                public int getVIPTrans() {
                        int aux=0;
                        for(int i=0; i<K;i++){
                                aux+=prop[i];
                        }
                                
                        return aux;
                }
                
                
                /**
                 * @return all OTHER transactions taking place
                 */
                public int getOTHERTrans() {
                        int aux=0;
                        for(int i=K; i<2*K;i++){
                                aux+=prop[i];        
                        }
                                
                        return aux;
                }
                 

                /**
                 * @param t is the type of transaction 
                 * @return the number of transactions taking place of a single type for VIP clients
                 */
                public int getVIPSingleTrans(int t) {
                        return prop[t];
                }

                
                
                /**
                 * @param t is the type of transaction
                 * @return the number of transactions taking place of a single type for OTHER users
                 */
                public int getOTHERSingleTrans(int t) {
                        return prop[K+t];
                }
                
                @Override
                public boolean isConsistent() {
                        return true;
                }

                @Override
                public String label() {
                        String stg = "";
                        stg+=" VIP transactions: ";
                                for(int k=0; k<K; k++){
                                            stg+=prop[k]+",";}
                        stg+=" OTHER Transactions: ";
                            for(int k=0; k<K; k++){
                                        stg+=prop[K+k]+",";}
                        
                        return stg + " VIPQ: " + getVIPQSize()+" OTHERQ: "+getOTHERQSize();
                }
                
        

        }
        /**
         * 
         * This class defines the events. 
         * An event has two components: Type which can have 4 values 
         * depending whether it represents an arrival of a VIP client, an 
         * arrival of OTHER user or a departure of: VIP client or OTHER user; and 
         * Transaction: defines which transaction is occurring.
         * 
         * @author Laura and Leonardo
         *
         */

        class BankEvent extends Event{
        
                static final int VIPARRIVAL = 1;
                static final int OTHERARRIVAL = 2; 
                static final int VIPDEPARTURE = 3;
                static final int OTHERDEPARTURE = 4;
                        int type; // ARRIVAL or DEPARTURE
                        int transaction;

        
                BankEvent(int type, int transaction) {
                this.type = type;
                this.transaction=transaction;}
        
        
                static EventsSet<BankEvent> getAllEvents(int K) {
                EventsSet<BankEvent> eSet = new EventsSet<BankEvent>();
        
                for (int i = 0; i < K; i++) {
                        eSet.add(new BankEvent(VIPARRIVAL, i));
                        eSet.add(new BankEvent(VIPDEPARTURE, i));
                        eSet.add(new BankEvent(OTHERARRIVAL, i));
                        eSet.add(new BankEvent(OTHERDEPARTURE,i));}
                
        
                return eSet;
        }

        /* (non-Javadoc)
         * @see java.lang.Object#toString()
         */
        @Override
                public String label() {
                String stg = "";
                switch (type) {
                        case (VIPARRIVAL) :
                                stg += "VIP arrival "+"transaction: "+(transaction+1);
                                break;
                        case (OTHERARRIVAL) :
                                stg += "OTHER Arrival "+"transaction: "+(transaction+1);
                                break;
                        case (VIPDEPARTURE) :
                                stg += "VIP Departure "+"transaction: "+(transaction+1);
                                break;
                        case (OTHERDEPARTURE) :
                                stg += "OTHER Departure "+"transaction: "+(transaction+1);
                                break;
                }
                return stg;
        }

}
        
        /**
         * This example describes a Bank office with a number of tellers N. There
         * is a queue for VIP clients (those with the bank's debit card), and
         * another one for OTHER users. 
         * 
         * @author Laura and Leonardo
         */
        class BankProcess extends SimpleMarkovProcess<BankState,BankEvent> {

                double[] lambda;//arrival rate for VIP (in 0) and OTHER(in 1)
                double[] mu;// service rates for transactions of type i  
                double[] Prob;//Probabilities for transaction types
                int M; // number of VIP servers 
                int S; // number of OTHER servers
                int K; // different transactions type
                int MaxVIPQ;//Max Capacity for VIP queue
                int MaxOTHERQ;//Max Capacity for OTHER queue
                private static final int VIPARRIVAL = BankEvent.VIPARRIVAL;
                private static final int OTHERARRIVAL = BankEvent.OTHERARRIVAL;
                private static final int VIPDEPARTURE = BankEvent.VIPDEPARTURE;
                private static final int OTHERDEPARTURE= BankEvent.OTHERDEPARTURE;
                
                
                /**
                 * @param lambda
                 * @param mu
                 * @param Prob
                 * @param M
                 * @param S
                 * @param K
                 * @param MaxVIPQ
                 * @param MaxOTHERQ
                 */
                public BankProcess(double[] lambda, double[] mu,double[] Prob, int M, int S, int K,int MaxVIPQ,int MaxOTHERQ) {
                        super(new BankState(M,S,K), BankEvent.getAllEvents(K));
                        this.M = M;
                        this.S = S;
                        this.K = K;
                        this.lambda = lambda;
                        this.mu = mu;
                        this.Prob=Prob;
                        this.MaxVIPQ=MaxVIPQ;
                        this.MaxOTHERQ=MaxOTHERQ;
                }
                /**
                 * Determines the active events. 
                 */
               
                @Override
                public boolean active(BankState i, BankEvent e) {
                        boolean result = false;
                        switch (e.type) {
                        case (VIPARRIVAL) ://VIPARRIVAL can occur if there is room in the VIP queue
                                                        {if(i.getVIPQSize()<MaxVIPQ){
                                                                        result =true;}
                                
                                break;
                        }
                         case (OTHERARRIVAL) ://OTHERARRIVAL can occur if there is room in the OTHER queue
                                        {if(i.getOTHERQSize()<MaxOTHERQ){
                                                result=true;}
                                break;
                         }
                         case (VIPDEPARTURE) ://if there is a VIP service in process of that transaction
                                        {if(i.getVIPSingleTrans(e.transaction)>0){
                                                result =true;}
                                break;
                         }
                         case (OTHERDEPARTURE) ://if there is an Other service in process of that transaction
                                        {if(i.getOTHERSingleTrans(e.transaction)>0){
                                                result =true;}
                        }
                        
                }
                        
                        return result;
        }

		/**
                 * Determines the possible destinations. 
                 */
                @Override
                public States<BankState> dests(BankState i, BankEvent e) {
                        StatesSet<BankState> set = new StatesSet<BankState>();
                        int[] props = new int[K+K+2];
                        for(int k=0;k<2*K;k++){
                        props[k] = i.getProp(k);}        //copy current values
                        
                        int VIPQ = i.getVIPQSize();
                        int OTHERQ = i.getOTHERQSize();
                        
                        switch (e.type) {
                        case (VIPARRIVAL) :
                                 if(i.getVIPTrans()>=M){// if there is no idle server, the client must stay in queue
                                        VIPQ++;
                                        set.add(new BankState(props,VIPQ,OTHERQ,K,M,S));} 
                                 else{
                                 props[e.transaction]++;
                                 set.add(new BankState(props,VIPQ,OTHERQ,K,M,S));} 
                                        
                        break;
                        case (OTHERARRIVAL):
                                 if(i.getOTHERTrans()>=S ){// if there is no idle server, the client must stay in queue
                                        OTHERQ++;
                                        set.add(new BankState(props,VIPQ,OTHERQ,K,M,S));} 
                                        
                                
                                        
                                 else{
                                                 props[K+e.transaction]++;
                                                 set.add(new BankState(props,VIPQ,OTHERQ,K,M,S));} 
                                        
                        break;
                                
                        case (VIPDEPARTURE) :
                                        if(i.getVIPQSize()==0 && i.getOTHERQSize()==0){
                                                        props[e.transaction]--;
                                                        set.add(new BankState(props,VIPQ,OTHERQ,K,M,S));}
                                                // if there is no VIPqueue and there is OTHERqueue, a VIP server can serve a normal client
                                        if(i.getVIPQSize()==0 && i.getOTHERQSize()>0){
                                                        props[e.transaction]--;
                                                        OTHERQ--;
                                                        for(int k=0; k<K; k++){
                                                                props[k]++;
                                                                set.add(new BankState(props,VIPQ,OTHERQ,K,M,S));        
                                                                props[k]--;}
                                                        }
                        
                                        if(i.getVIPQSize()>0){
                                                        props[e.transaction]--;
                                                                VIPQ--;
                                        for(int k=0; k<K; k++){
                                                        props[k]++;
                                                        set.add(new BankState(props,VIPQ,OTHERQ,K,M,S));        
                                                        props[k]--;}
                        }
                        break;
                        case (OTHERDEPARTURE) :
                                if(i.getOTHERQSize()==0){
                                        props[K+e.transaction]--;
                                        set.add(new BankState(props,VIPQ,OTHERQ,K,M,S));}
                                
                                if(i.getOTHERQSize()>0){
                                        props[K+e.transaction]--;
                                        OTHERQ--;
                                for(int k=0; k<K; k++){
                                        props[K+k]++;
                                        set.add(new BankState(props,VIPQ,OTHERQ,K,M,S));
                                        props[K+k]--;}
                                }
                        break;        
                        
            }
                        
                        return set;
                        

                }

		/**
                 * Determines the rates. 
                 */
                @Override
                public double rate(BankState i, BankState j, BankEvent e) {
                        double result = 0;
                        
                        switch (e.type) {
                                        case (VIPARRIVAL) :
                                                result = lambda[0]*Prob[e.transaction];  
                                        
                                        break;
                                        
                                        case (OTHERARRIVAL) :
                                                result = lambda[1]*Prob[e.transaction];  
                                        break;
                                                
                                        case (VIPDEPARTURE) :
                                                result=mu[e.transaction];
                                        break;

                                        case (OTHERDEPARTURE) :
                                                result=mu[e.transaction];
                                        break;
                                                
                                                         }
                        return result;
                                                }
                                        
                        
                

                @Override
                public String description() {
                        
                        String stg = "Bank with VIP Servers SYSTEM\n\n";
                        stg += "Multiple server queue with " + this.M + " servers for VIP clients\n";
                        stg+="and "+this.S+" servers for normal clients\n\n";
                        try {
                                                        stg+="VIP clients mean waiting time: "+((this.getMOPsAvg("VIPQueue Length")/this.lambda[0])*60)+" minutes\n\n\n";
                                                } catch (NotUnichainException e) {
                                                        e.printStackTrace();
                                                }
                        if(this.S==0){
                                stg+="Time desired cannot be achieved with that number of servers!!!!\n";
                                stg+="Try adding more servers";
                        }
                        
                        return stg;
                }

        }

                class BankVIP{ 

                /**
                 * @param a
                 * @throws NotUnichainException
                 */
                public static void main(String[] a) throws NotUnichainException {
                double[] lambda = new double[2];
                double[] prob = new double[0];
                double[] mu=new double[0];
                int M=0;
                int S=0;
                int N=0; 
                int K=0;
                int MaxVIPQ=0;
                int MaxOTHERQ=0;
                double t=0;
                double T=0;
                
                
               
                    BufferedReader rdr =
                                new BufferedReader(new InputStreamReader(System.in));
                        try {
                                System.out.println("VIP Arrival Rate: ");
                                lambda[0] = Double.parseDouble(rdr.readLine());
                                System.out.println("Other Arrival Rate: ");
                                lambda[1] = Double.parseDouble(rdr.readLine());
                                System.out.println("Total Serves: ");
                                N = Integer.parseInt(rdr.readLine());
                                System.out.println("Total Transaction types: ");
                                K = Integer.parseInt(rdr.readLine());
                                mu=new double[K];
                                for(int k=0;k<K;k++ ){
                                        System.out.println("Service Rate for transaction "+(k+1)+":");
                                        mu[k]=Double.parseDouble(rdr.readLine());}
                                prob= new double[K];
                                for(int k=0;k<K;k++ ){
                                        System.out.println("Probability for transaction "+(k+1)+":");
                                        prob[k]=Double.parseDouble(rdr.readLine());}
                                
                                System.out.println("Queue Capacity: ");
                                MaxVIPQ = Integer.parseInt(rdr.readLine());
                                System.out.println("Desired waiting time for VIP users in minutes:");
                                    T=Double.parseDouble(rdr.readLine());
                                    T=T/60;
                                    
                                
                                
                        } catch (IOException e) {
                        };
                
                
                
                
                t=T+1;
                S=N;
                MaxOTHERQ=MaxVIPQ;
                
                for (M=0;M<N;M++){
                       if(t>T){
                    BankProcess TheModel = new BankProcess(lambda,mu, prob,  M,  S,  K,MaxVIPQ, MaxOTHERQ );
                
                    TheModel.generate();
                t=TheModel.getMOPsAvg("VIPQueue Length")/lambda[0];
                S--;

                
                                       }
                       
                       else{break;}
                   }
               
               BankProcess TheModel = new BankProcess(lambda,  mu, prob,  M,  S,  K,MaxVIPQ, MaxOTHERQ );
               TheModel.printAll();
               TheModel.showGUI();
                                             
               

        }

}




