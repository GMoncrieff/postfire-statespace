## A state-space model of vegetation activity and fire recovery in the Cape Floristic Region

The model uses the folllowing state update equation:

$$z_{i,t} = (z_{i,t-1} + z_{i,t-1}\lambda_i(1-\frac{z_{i,t-1}}{\gamma_i})) (1-q_{i,t-1}) + q_{i,t-1}\alpha_i +  A_i\sin(x_t+\phi_i) +  \epsilon_{i,t}$$
$$\epsilon_{i,t} \sim \mathcal{N}(0,\sigma_p^{2})$$
$$x_{i,t} \sim \mathcal{bernoulli}(q_{i,t})$$
$$y_{i,t} \sim \mathcal{N}(z_{i,t},\sigma_o^{2})$$

This represents logistic growth with seasonality and a fire switch that resets the vegetation state when fire occurs.
$z_{i,t}$ is the hidden state for site $i$ at time $t$
 $x_{i,t}$ is the observed occurence of fire and is observed with error, hence it is sampled from a bernoulli distribution wih a rate parameter $q_{i,t}$
 $A_i\sin(x_t+\phi_i)$ introduces seasonality. 
 $\lambda, \gamma, \alpha, A$ and $\phi$ are the parmaeters for growth rate, max ndvi, ndvi immediately after fire, seasonal amplitude, and seasonal phase shift respectively

An example time series below shows the median estiatmed state (black line) with 90% CI (dark shaded area) and observations (red points). Forecasted future states are shown in the light gray shaded area

![example ts](test_ts.png)