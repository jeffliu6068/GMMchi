def solve_gaussians(m1, m2, cov1, cov2, s1, s2, minimum, maximum):
    std1 = np.sqrt(cov1)
    std2 = np.sqrt(cov2)

    a = 1/(2*std1**2) - 1/(2*std2**2)
    b = m2/(std2**2) - m1/(std1**2)
    c = m1**2 / (2*std1**2) - m2**2 / (2*std2**2) - np.log((std2*s1)/(std1*s2))
    # GET THE INTERSECTION THAT IS BETWEEN THE TWO MEANS
    return np.sort([x for x in np.roots([a[0], b[0], c[0]]) if x < maximum and x > minimum])
