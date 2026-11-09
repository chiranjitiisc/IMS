import numpy as np
rn = np.random.default_rng()

# rolling a die
# method_1
x  = np.array([1, 2, 3, 4, 5, 6])                    # every possible value
px = np.array([0.05, 0.05, 0.10, 0.10, 0.20, 0.50])  # Pr(x_i), adds up to 1
gx = x**2                                            # g(x_i)

total = 0.0
for i in range(len(x)):      # visit each value exactly once
    total += gx[i] * px[i]   # g(x_i) * Pr(x_i)
print(total)                 # 25.75, the same every run

# method 2
M = 20
samples = rn.choice(x, size=M, p=px)    # faces drawn with the loaded probabilities
print(np.mean(samples**2))              # about 25.75; gave 28.25 for M = 20
print(samples)