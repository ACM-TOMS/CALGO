/*
 * Created on 21/10/2005
 */
package jmarkov.basic;

/**
 * This class represents the joint information of a value function and a policy
 * which summarizes the solution to a problem.
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 *
 * @param <S> state
 * @param <A> action
 */
public class Solution <S extends State, A extends Action>{
    private ValueFunction<S> valueFunction;
    private Policy<S,A> policy;

    /**
     * Builds a solution given a value funtcion and a policy
     * @param valueFunction value function
     * @param policy policy
     */
    public Solution(ValueFunction<S> valueFunction, Policy<S,A> policy) {
        super();
        this.valueFunction = valueFunction;
        this.policy = policy;
    }

    /**
     * Returns the Policy.
     * @return Returns the policy.
     */
    public Policy<S,A> getPolicy() {
        return policy;
    }

    /**
     * Returns the valueFunction.
     * @return Returns the valueFunction.
     */
    public ValueFunction<S> getValueFunction() {
        return valueFunction;
    }

}
