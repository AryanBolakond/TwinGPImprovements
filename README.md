Although TwinGP provides an efficient framework for large-scale GP modeling, its performance depends on a number of design choices.

Many of these decisions are currently guided by intuition, prior experience,
or trial and error. For example, global points may be selected using twinning,
random sampling, or clustering-based methods. The number of global and local
points can vary substantially, affecting both prediction accuracy and computational cost. Similarly, different kernel configurations may lead to different
predictive behavior across datasets.
Despite the success of TwinGP, there is currently no systematic framework
for identifying which combinations of design choices produce the best tradeoff
between prediction accuracy and computational efficiency. This motivates the
central question of this project: Which TwinGP design choices lead to the most
accurate predictions while minimizing computational cost?

The resulting experiments will allow us to quantify the effect of individual
design choices on prediction accuracy and computational cost. Based on these
results, we will develop recommendations for configuring TwinGP to achieve
strong predictive performance with minimal runtime.
