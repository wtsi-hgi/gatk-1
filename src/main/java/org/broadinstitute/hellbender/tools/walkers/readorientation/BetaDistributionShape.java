package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.utils.Utils;

public class BetaDistributionShape {
    private double alpha;
    private double beta;

    public BetaDistributionShape(final double alpha, final double beta){
        Utils.validateArg(alpha > 0 , "alpha must be greater than 0 but got " + beta);
        Utils.validateArg(beta > 0, "beta must be greater than 0 but got " + beta);

        this.alpha = alpha;
        this.beta = beta;
    }

    public double getAlpha() {
        return alpha;
    }

    public double getBeta() {
        return beta;
    }
}
