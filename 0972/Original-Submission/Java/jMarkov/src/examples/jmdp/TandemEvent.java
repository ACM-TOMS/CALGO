package examples.jmdp;

import jmarkov.basic.Event;

/**
* This class represents an event in a tandem queuing system it is used by the Access Control examples
   * @author Daniel Silva 
     * 
     */
  

public class TandemEvent extends Event {
   public enum TEvent {
        /** When a customer arrives */
        Arrival,
        /** When a customer completes service at station 1 */
        Service1,
        /** When a customer completes service at station 2 */
        Service2,
    }

    private TEvent e;

    /**
     * @param e
     *            An event of the types listed
     */
    public TandemEvent(TEvent e) {
        this.e = e;
    }

    /**
     * @return An event of the types listed
     */
    public TEvent get() {
        return e;
    }

    @Override
    public String label() {
        return e.name();
    }

    @Override
    public String description() {
        return "Event : " + e.name();
    }

    @Override
    public int compareTo(Event e1) {
        if (!(e1 instanceof TandemEvent))
            return -1;
        TandemEvent e2 = (TandemEvent) e1;
        return e.compareTo(e2.e);
    }
}