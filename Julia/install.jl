using Pkg

# Update Package State
Pkg.update()

# Necessary external modules
Pkg.add("LsqFit")        # Allows Least Square Fits
Pkg.add("FFTW")          # Makes Fast Fourrier Tranforms
Pkg.add("Statistics")    # contains statistics functions
Pkg.add("Bootstrap")     # Allows Bootstrap

# Update Package State
Pkg.update()
