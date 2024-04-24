# Reactivity-Ratios
Software to calculate reactivity ratios in radical copolymerization with the use of the terminal model

There are several codes that can perform the task a) Contour (freeware), a stand-alone program which alows for the IUPAC recommended method, but also with additional features like fitting propagation rate coefficients with the penultimate model and fitting the Arrhenius equation b) an Excel sheet where the IUPAC recommended method and c) open source Python code doing the same. All three applications use the same file format; f0, X, F, ∆F in a comma separated (.CSV) file.
Contour can be dowloaded via this link [Contour](https://drive.google.com/file/d/1lYicEcPbuQiohQGE4x5RUVALTVU9sWnL/view?usp=drive_link) . The other applications are here in the repository. In the Contour folder you can also find instruction videos for Contour and for the Excel file converter from f to F.
In the IUPAC recommended method the f0-X-F plane is analysed. The experiments are set-up by selecting an initial feed fraction of the monomer (f0) and performing a copolymerization experiment. The cumulative copolymer composition (F) is analysed and the total molar monomer conversion (X, for both monomers together and values between 0 and 1) at the moment of copolymer isolation is determined. An error in F needs to be estimated (∆F) and entered into the calculation, this value can used as a weighing factor in the fit.  Alternatively one can use uniform weights of all data (choice of weighting scheme).
The IUPAC working group advises experimentalists to run several experiments at different f0 values and always determine the conversion. The data of f0-F (low conversion) experiments can still be used here, assuming that the conversion is also measured and present in the .CSV file. 
With the IUPAC supported approach, datapoints on the f0-F-X plane are used to determine the reactivity ratios and, as long as the conversion is known, low and high conversion experiments can be used and mixed. 
Regarding higher conversion data, one has to be aware that sometimes one of the two monomers is completely consumed. These data could be taken out of the set because now homopolymer is produced only and no new information is added to the dataset. 

Which error scheme is to be used?
If all the data have equal weight, uniformous weighing should be used, but in general it is preferred to use weighing by the error in each individual F.
If the error in X is substantial and we are not dealing with F-values calculated from f with proper error propagation (see below), the error in all variables method should be used (EVM, currently not covered).
In some case one has several experiments at different f0 where one is following the feed composition with conversion (f-X). Those data can be converted to f0-X-F data with a separate Excel sheet. The errors in f (and f0) and X each are propagated into an error in F which needs to be taken into account in the fit as a weighing factor (individual error per point). As the errors in X and f are already error propagated in an error in F, EVM should not be used.
The program uses numerical integration of the copolymerization equation and visualization of the sum of squares of residuals to find the optimal values of the reactivity ratios. Each measurement is weighed in the fit with the use of the error in F (∆F). It can take minutes to display the solution, due to the numerical approach. One could start with a relatively large search range for r1 and r2 and then zoom in nearer to the minimum. In that way the joint confidence interval and the r-values get more accurate. In Contour this is partially automatic, still it is efficient to already define a search range by hand. The errors are often not symmetric around the optimum because with our visualization approach we determine the true shape of the joint confidence interval and not the ellipsoid  approximation. 
