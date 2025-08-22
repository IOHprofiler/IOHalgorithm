import numpy as np
import ioh
import math
# ------------------------ (1+1)EA ------------------------

class OnePlusLambdaLamdbdaEA:
    lam: int = 10
    budget_factor = 5
    b = 2.0/3.0
    a = pow(1.5,0.25)
    seed = 42
    
    rng = None

    def sample_conditional_binomial(self, p: float, n: int) -> int:
        """Sample Binomial(n, p), but ensure in [1, n]."""
        k = sum(1 for _ in range(n) if self.rng.random() < p)
        if k < 1:
            k = 1
        elif k > n:
            k = n
        return
    
    def __call__(self, f, seed = 1) -> tuple[float, np.ndarray]:
        """
        Run the (1+λ)EA on a single IOH problem object `f`.

        Returns
        -------
        (fx_best, x_best) : tuple[float, np.ndarray]
            Best fitness found and corresponding solution.
        """
        self.seed = seed
        self.rng = np.random.default_rng(self.seed)

        n = f.meta_data.n_variables
        budget = self.budget_factor * n * n

        # Initialize x within bounds (PBO: binary {0,1})
       
        x = self.rng.integers(0,2,n)
        fx = f(x)


        while f.state.evaluations < budget:
            # Generate λ offspring
            best_m_y = None
            best_m_fy = None


            update_lambda_flag = False
            
            mutation_rate = self.lam / n
            crossover_rate = 1 / self.lam
            ell = self.sample_conditional_binomial(mutation_rate,n)

            for _ in range(int(self.lam)):
                
                y = x.copy()
                idx = self.rng.choice(n, size=ell, replace=False)
                # bit flip for {0,1}: new = lb + ub - old (works for 0/1 bounds)
                y[idx] = 1 - y[idx]
                fy = f(y)

                if best_m_y is None:
                    best_m_y, best_m_fy = y, fy
                else:
                    if fy > best_m_fy:
                        best_m_y, best_m_fy = y, fy


            for _ in range(int(self.lam)):
                update_flag = False
                y = x.copy()
                for i in range(n):
                    if self.rng.uniform() < crossover_rate:
                        y[i] = best_m_y[i]
                        if x[i] != best_m_y[i]:
                            update_flag = True
                
                if not update_flag:
                    fy = fx
                elif np.all(y == best_m_y):
                    fy = best_m_fy
                else:
                    fy = f(y)


                if fy > fx:
                    update_lambda_flag = True
                if fy >= fx:
                    x, fx = y, fy
            
            if best_m_fy > fx:
                update_lambda_flag = True
                x, fx = best_m_y, best_m_fy
                
            if fx >= f.optimum.y:
                break

            if update_lambda_flag:
                self.lam = max(self.lam * self.b, 1)
            else:
                self.lam = min(self.lam * 2, n)
        return fx, x