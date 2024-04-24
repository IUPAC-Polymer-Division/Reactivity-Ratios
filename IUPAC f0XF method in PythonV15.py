# version 04-04-2024, Disclaimer, no guarantee is given that this code works correcty
import sys
# Global variables
f10Step=0 # for optimizing f10 if needed
f10Start=0 # "
f10End=0 # "
f100=0 #"
f100min=0 #"
f10old=0 #""
sdr1=0 # standard deviation in r1
sdr2=0 # standard deviation in r2
l90 = 0 # 90% confidence level for JCI
l95 = 0 # 90% confidence level for JCI
rtheor = 0 # theoretical ssr based on errors from user
stepsizeconv = 0 # stepsize for numerical intergration of copo eq. Now 0.001, 0.003 also works faster, less accurate, we use 0.003 when optimizing f10
r1range = 0 # search range r1 (1/2 this range used around the centre of search region)
r2range = 0 # search range r2 (1/2 this range used around the centre of search region)
nosample = 0 # number of measurements
convs = 0 # conversion in integration loop
fnew = 0 # updated small f1 in copo eq
f2new = 0 # updated small f2 in copo eq
imaxi = 0 # maximum number of integration steps
stepsizer1 = 0 #stepsize for r1, determines the search grid
stepsizer2 = 0 #stepsize for r2, determines the search grid
k = 0 # counter for SSR array x direction
j = 0 # counter for SSR array y direction
delta = 0 # residual for one measurement
errscheme = 0 #errscheme=1 using errors given by user, errscheme=2 all measurements weight of 1
r1or = 0 # centre of search region for r1
r2or = 0 # centre of search region for r2
sample = 0 # counter for measurements
answer = 0 # used in question
i=0 # counter
path = ""
dummy = ""
XYnaam = ""
F10 = 0 # big F1 at start of the reaction
F1avest = 0 # estimated big F1 in one intergration step
F1totest = 0 # is the smae as Fcum, total big F after integration over conversion (has to be divided by total no of integration steps imaxi)
r1 = 0 #r1
r2 = 0 #r2
convber = 0 # calculated conversion
deltatotal = 0 # total SSR
deltamin = 0 # lowest value in SSR array
r1min = 0 # corresponding optimal value for r1
r2min = 0 # corresponding optimal value for r2

# Arrays
f10 = [0] * 201 # f10
F = [0] * 201 # big F as measured
deltaF = [0] * 201 #error in big F
conv = [0] * 201 #conversion
f1 = [0] * 201 #small f, monomer ratio in monomer mixture
w = [0] * 201 # weight of meaurement
remark = [""] * 201 # remark about whether this datapoint has a reasonable error or bigger than error estimate or one of the monomers run out
F1tot = [0] * 201 # array for big F after integration
Fopt = [0] * 201 # array of big F at optimum values of r1 r2

import numpy as np
ssr = np.empty((102, 102), dtype=float) # sum of squares of residuals array, used to find JCI, can also be saved
def Numint():# numerical integration of the copolymerization equation. with the same f10 one integration with stops at the conversions
    global deltatotal, F1tot
    deltatotal=0
    f10old=100
    for sample in range(1, nosample + 1):
   
        fnew = f10[sample]
        f2new = 1 - fnew
        if f10old != f10[sample]:
            f10old=f10[sample]
            convs = 0
            F1totest = (r1 * fnew ** 2 + fnew * f2new) / (r1 * fnew ** 2 + 2 * fnew * f2new + r2 * f2new ** 2) # first value has already to be calculated
            F1avest = 0
            imaxi = round(conv[sample] / stepsizeconv)  # total number of integration steps
            i_val = 0
        else: # if we are within the same f10 series we continue to integrate for the differenct conversions, we already sorted according to f10 and conv
            i_val=imaxi
            imaxi = round(conv[sample] / stepsizeconv)  # total number of integration steps for this conversion
            
        while i_val != imaxi:  # numerical integration
            i_val += 1
            convs += stepsizeconv
            fnew = (f10[sample] - (F1totest / i_val)*convs) / (1 - convs)
            f2new = 1 - fnew
            F1totest += F1avest # summation over all conversion steps
            F1avest = (r1 * fnew ** 2 + fnew * f2new) / (r1 * fnew ** 2 + 2 * fnew * f2new + r2 * f2new ** 2) # F over one integration step

            if i_val == imaxi:
                F1tot[sample] = F1totest / imaxi  # F1tot is the result, as we took imaxi steps we need to divide by imax

        delta = F1tot[sample] - F[sample]
        deltatotal += w[sample] * delta ** 2
        
        
import tkinter as tk
from tkinter import filedialog

def open_file_dialog():
    global remark, f10, conv, F, deltaF, f1, nosample  # Declare variables as global
    root = tk.Tk()
    root.withdraw()  # Hide the main window

    file_path = filedialog.askopenfilename(title="Select File with f10,X,F, delta F data")
    if not file_path:  # User pressed cancel or closed the dialog
        print("File selection cancelled.")
        return 0  # 
      
    if not file_path.lower().endswith('.csv'):
        print("Selected file is not a CSV file.")
        return 0  # Or handle the error accordingly

    if file_path:
        print("Selected file:", file_path)
        with open(file_path,'r') as file:
            i = 0
            while True:
                i += 1
                line = file.readline()
                if not line:
                    break
            
                data = line.split(',')
                remark[i] = ""
                f10[i] = float(data[0])
                conv[i] = float(data[1])
                if conv[i]<=0.003:
                    conv[i]=0.003 # needed to assure at least one loop trough the integration, 0.003 is the highest stepsize at the moment (single f10 set)
                    
                F[i] = float(data[2])
                deltaF[i] = float(data[3])
                f1[i] = (f10[i] - F[i]*conv[i]) / (1 - conv[i])
                if f1[i] < 0:
                    f1[i] = 0.000
                    remark[i] = "Seems no more monomer 1, consider removing,"
                    
                if f1[i] >= 1:
                    f1[i] = 1.000
                    remark[i] = "Seems no more monomer 2, consider removing,"

        return i
def sort_arrays(f10, conv, F, deltaF, remark, f1):
        # Perform bubble sort based on f10 and conv
    for i in range(1, nosample- 1):
        for j in range(i + 1, nosample):
            if f10[i] > f10[j] or (f10[i] == f10[j] and conv[i] > conv[j]):
                # Swap elements
                f10[i], f10[j] = f10[j], f10[i]
                conv[i], conv[j] = conv[j], conv[i]
                F[i], F[j] = F[j], F[i]
                deltaF[i], deltaF[j] = deltaF[j], deltaF[i]
                remark[i], remark[j] = remark[j], remark[i]
                f1[i], f1[j] = f1[j], f1[i]
                # Using tuple unpacking for swapping multiple values in one line


def get_inputs():
    global errscheme, r1or, r2or, r1range, r2range, nosample  # Declare variables as global
    global save_jci, save_residuals, save_ssr  # Declare global boolean variables
    
    # Retrieve input values from entry fields
    errscheme = int(entry1.get())
    r1or = float(entry2.get())
    r2or = float(entry3.get())
    r1range = float(entry4.get()) / 2.0
    r2range = float(entry5.get()) / 2.0

    # Display the inputs
    print("Errorscheme:", errscheme)
    print("r1 estimate:", r1or)
    print("r2 estimate:", r2or)
    print("r1range: +-", r1range)
    print("r2range: +-", r2range)

    # Update global boolean variables based on user's choices
    save_jci = True if jci_choice.get() == "Yes" else False
    save_residuals = True if residuals_choice.get() == "Yes" else False
    save_ssr = True if ssr_choice.get() == "Yes" else False
      
    # Close the pop-up window
    popup_window.destroy()

    # Open file dialog
    nosample = open_file_dialog() - 1 # number of measurements
    if nosample==-1:
        sys.exit()
        
    print(f"{nosample} datapoints")

    # Close the root window to terminate the main loop
    root.quit()

# Create a Tkinter window

root = tk.Tk()
root.withdraw()  # Hide the main window

# Create a pop-up window
popup_window = tk.Toplevel(root)
popup_window.title("Input Parameters")

# Create labels and entry fields
label1 = tk.Label(popup_window, text="Weighting based on given errors in F (1) or uniform weighting (2)?")
label1.grid(row=0, column=0, padx=5, pady=5)
entry1 = tk.Entry(popup_window)
entry1.grid(row=0, column=1, padx=5, pady=5)

label2 = tk.Label(popup_window, text="Middle of search region for r1:")
label2.grid(row=1, column=0, padx=5, pady=5)
entry2 = tk.Entry(popup_window)
entry2.grid(row=1, column=1, padx=5, pady=5)

label3 = tk.Label(popup_window, text="Middle of search region for r2:")
label3.grid(row=2, column=0, padx=5, pady=5)
entry3 = tk.Entry(popup_window)
entry3.grid(row=2, column=1, padx=5, pady=5)

label4 = tk.Label(popup_window, text="Search range around r1:")
label4.grid(row=3, column=0, padx=5, pady=5)
entry4 = tk.Entry(popup_window)
entry4.grid(row=3, column=1, padx=5, pady=5)

label5 = tk.Label(popup_window, text="Search range around r2:")
label5.grid(row=4, column=0, padx=5, pady=5)
entry5 = tk.Entry(popup_window)
entry5.grid(row=4, column=1, padx=5, pady=5)

# Label for file saving options
save_label = tk.Label(popup_window, text="Save Files")
save_label.grid(row=5, columnspan=3, padx=5, pady=5)

# Dropdown menus for file saving options
jci_label = tk.Label(popup_window, text="JCI:")
jci_label.grid(row=6, column=0, padx=5, pady=5)
jci_choice = tk.StringVar(popup_window)
jci_choice.set("No")
jci_option_menu = tk.OptionMenu(popup_window, jci_choice, "Yes", "No")
jci_option_menu.grid(row=6, column=1, padx=5, pady=5)

residuals_label = tk.Label(popup_window, text="Residuals Plot:")
residuals_label.grid(row=7, column=0, padx=5, pady=5)
residuals_choice = tk.StringVar(popup_window)
residuals_choice.set("No")
residuals_option_menu = tk.OptionMenu(popup_window, residuals_choice, "Yes", "No")
residuals_option_menu.grid(row=7, column=1, padx=5, pady=5)

ssr_label = tk.Label(popup_window, text="SSR Surface:")
ssr_label.grid(row=8, column=0, padx=5, pady=5)
ssr_choice = tk.StringVar(popup_window)
ssr_choice.set("No")
ssr_option_menu = tk.OptionMenu(popup_window, ssr_choice, "Yes", "No")
ssr_option_menu.grid(row=8, column=1, padx=5, pady=5)

# Submit Button
submit_button = tk.Button(popup_window, text="Load File", command=get_inputs)
submit_button.grid(row=9, columnspan=3, padx=5, pady=10)
# Set initial values for input parameters
entry1.insert(0, "1")  # Initial value for weighting scheme
entry2.insert(0, "0.4")  # Initial value for r1 estimate
entry3.insert(0, "0.6")  # Initial value for r2 estimate
entry4.insert(0, "0.4")  # Initial value for search range around r1
entry5.insert(0, "0.4")  # Initial value for search range around r2

jci_choice.set("No")  # Initial value for JCI dropdown menu
residuals_choice.set("No")  # Initial value for Residuals Plot dropdown menu
ssr_choice.set("No")  # Initial value for SSR Surface dropdown menu

# Run the Tkinter event loop
root.mainloop()
sort_arrays(f10, conv, F, deltaF, remark, f1) # sorting according to f10 and conversion for nicer table and numeric integration
r1 = r1or
r2 = r2or
# Fitting reactivity ratios based on differences in F

'filling up the weights and calculating theor SSR'
rtheor = 0
for sample in range(1, nosample + 1):
    if errscheme == 1:
        w[sample] = 1 / deltaF[sample] ** 2  # simple weighing according to deltaF
        rtheor += w[sample] * deltaF[sample] ** 2
    else:
        w[sample] = 1  # uniform weighing
        rtheor += w[sample] * deltaF[sample] ** 2

f10loop = False
AreAllElementsSame = True
for i in range(2, nosample + 1):
    if f10[i] != f10[1]:
        AreAllElementsSame = False

response = None
if AreAllElementsSame:
    response = input("One set of f10. Do you want to optimize f10 (+- 0.005) ? (calc. longer) [Yes/No]: ")
    response = response.lower()

if response == "yes":
    f10loop = True
    f10Step = 0.001
    f10Start = f10[1] - 0.005
    f10End = f10[1] + 0.005
 
   

print()
print("Calculation in progress, wait a few minutes.......")
print()
if f10loop:
 print(f"Looking for an optimum in f10 between: {round (f10Start, 4)} and {round(f10End,4)}")
print(f"Looking for an optimum of r1 between: {round(r1or - r1range,4)} and {round(r1or + r1range,4)}")
print(f"Looking for an optimum of r2 between: {round(r2or - r2range,4)} and {round(r2or + r2range,4)}")
r1min = 10000
r2min = 10000
deltamin = 10000
stepsizeconv = 0.001  # stepsize for integration in the actual calculation, combined with the fineness of the grid (now 100 X 100) this determines the calculation time
stepsizer1 = r1range * 2.0000 / 100.0000  # fineness of the grid we are looking for r1
stepsizer2 = r2range * 2.0000 / 100.0000  # fineness of the grid we are looking for r2
if not f10loop:
 k = -1  # counter for SSR array r1 axis
 r1 = r1or - r1range - stepsizer1
 while not r1 > r1or + r1range:
    r1 += stepsizer1
    k += 1  
    r2 = r2or - r2range - stepsizer2
    j=-1 # counter for SSR array r2 axis
    while not r2 > r2or + r2range:
        r2 += stepsizer2
        j += 1
             
        Numint()

        ssr[k,j] = deltatotal # filling SSR array, for some reason x and y has to be switched in making the contour line later
        
        if deltatotal < deltamin:  # remember the lowest values
            deltamin = deltatotal
            r1min = r1
            r2min = r2
            Fopt = F1tot.copy()  # Copy the values of F1tot to Fopt
else:   # extra loop for f10   
    stepsizeconv = 0.003  # stepsize for integration in the actual calculation for f10 optimization, slightly less accurate to save time, watch out for conv<stepsizeconv
    f100min = 1  
    f100 = f10Start
    f10old=f10[1]
    for i in range(1 ,nosample+1):
        f10[i]=f100
    
    while f100 <= f10End: 
       
        r1 = r1or - r1range - stepsizer1
        while r1 <= r1or + r1range:  
            r1 += stepsizer1
          
            r2 = r2or - r2range - stepsizer2
           
            while r2 <= r2or + r2range:
                r2 += stepsizer2
                                           
                Numint()
                              
                if deltatotal < deltamin:  # remember the lowest values
                    deltamin = deltatotal
                    r1min = r1
                    r2min = r2
                    Fopt = F1tot.copy()  # Copy the values of F1tot to Fopt
                    f100min = f100

                               
        # Increment f100 in the outermost loop
        f100 += f10Step
print("\n")
print("Calculation finished")      

if f10loop:
    r1min = 10000
    r2min = 10000
    deltamin = 10000
    print()
    print("Your optimal f10 value is:", f100min , "(was ", f10old, ")")  
    if f100min == f10Start or f100min == f10End:
        print("Your optimal f10 value is on the edge of search region!:", round(f100min, 4))  
    for i in range(1, nosample + 1):
        f10[i] = f100min
        f1[i] = (f10[i] - F[i] * conv[i]) / (1 - conv[i]) # need to recheck the f1 values with the new f10
        if f1[i] < 0:
                    f1[i] = 0.000
                    remark[i] = "Seems no more monomer 1, consider removing,"
                    
        if f1[i] >= 1:
                    f1[i] = 1.000
                    remark[i] = "Seems no more monomer 2, consider removing,"
                    
      # filling SSR array after establishing optimal f10  
    stepsizeconv = 0.003 # stepsize does not need to be the same for f10 search loop and final calculation, minimum is obtained again anyway
    k = -1  # counter for SSR array r1 axis
    r1 = r1or - r1range - stepsizer1
    while r1 <= r1or + r1range:  
              r1 += stepsizer1
              k += 1  
              r2 = r2or - r2range - stepsizer2
              j = -1  # counter for SSR array r2 axis
              while r2 <= r2or + r2range:
                  r2 += stepsizer2
                  j += 1
                              
                  Numint()


                  ssr[k,j] = deltatotal # filling SSR array after establishing optimal f10   
                  if deltatotal < deltamin:  # remember the lowest values, in case stepsizeconv in the main loop is different from this loop we need to recalculate the minima
                      deltamin = deltatotal
                      r1min = r1
                      r2min = r2
                      Fopt = F1tot.copy()  # Copy the values of F1tot to Fopt
                      
          
        

r1 = r1min
r2 = r2min
if f10loop:
    l95 = deltamin + 7.815 * rtheor / (nosample - 3) # using Chi squared distribution with 3 degrees of freedom
    l90 = deltamin + 6.251 * rtheor / (nosample - 3) # not used now, can be set as alternative level for JCI contourline later
else:
     l95 = deltamin + 5.99147 * rtheor / (nosample - 2) # using Chi squared distribution with 2 degrees of freedom
     l90 = deltamin + 4.60517 * rtheor / (nosample - 2) # not used now, can be set as alternative level for JCI contourline later

if r1 <= r1or - r1range +stepsizer1 or r1 >= r1or + r1range-stepsizer1:
    print(f"Optimum r1 ({r1}) on the edge of search range, shift search range!, restart")
    #stop the program
    sys.exit()
if r2 <= r2or - r2range+stepsizer2 or r2 >= r2or + r2range-stepsizer2:
    print(f"Optimum r2 ({r2}) on the edge of search range, shift search range!, restart")
    #stop the program
    sys.exit()
print()
print ("JCI based on Chi squared distribution combined with the user given errors in F")
print(f"95% confidence limit at contour height {round(l95, 4)}  90% confidence limit at contour height {round(l90, 4)}")
print()

# show JCI and minimum
import numpy as np
import matplotlib.pyplot as plt


# Transpose ssr to have r1 along the x-axis !!
ssr = np.transpose(ssr)
# crop array to take away noisy edges, else spoiling graph and the determination of outer limits of JCI to determine Standard deviations, if needed increase to 2%
crop_percentage=1
# Calculate cropping indices
crop_start = int(ssr.shape[0] * crop_percentage / 100)
crop_end = int(ssr.shape[0] * (100 - crop_percentage) / 100)

# Crop ssr array
ssr_cropped = ssr[crop_start:crop_end, crop_start:crop_end]

# Create a 2D surface based on the cropped data
x = np.linspace(r1or - r1range, r1or + r1range, ssr_cropped.shape[0])
y = np.linspace(r2or - r2range, r2or + r2range, ssr_cropped.shape[1])
X, Y = np.meshgrid(x, y)

# Specify the desired contour value
contour_value = l95

# Plot the 2D surface with only one contour line at the specified value
contour = plt.contour(X, Y, ssr_cropped, levels=[contour_value], colors='k')
plt.xlabel('$r_{1}$')
plt.ylabel('$r_{2}$')
plt.title('Joint Confidence Interval at 95% probability')

# Extract contour vertices
# Check if contour collections exist before accessing them
if contour.collections:
    contour_vertices = contour.collections[0].get_paths()[0].vertices
else:
    print("No contour found in search region")

# Find highest and lowest x values
highest_x, lowest_x = np.max(contour_vertices[:, 0]), np.min(contour_vertices[:, 0])
# Find highest and lowest y values
highest_y, lowest_y = np.max(contour_vertices[:, 1]), np.min(contour_vertices[:, 1])

# Plot the minimum point
plt.scatter(r1, r2, color='red', marker='o', label='Optimum')

# If JCI is OK, calculate standard deviations from top outer limits of contourline, just rough estimate, rather use JCI itself
if highest_x<r1or+r1range and lowest_x>r1or-r1range and highest_y<r2or+r2range and lowest_y>r2or-r2range:
    sdr1 = (highest_x - lowest_x) / 4.9
    sdr2 = (highest_y - lowest_y) / 4.9
    print("Maybe zoom in further")
else:
    print("JCI not fully in range, SDs set to zero")
    sdr1 = 0
    sdr2 = 0
    
if sdr1 < 0 or sdr2 < 0:
    print("No reliable JCI found in your search range, SDs set to zero")
    sdr1 = 0
    sdr2 = 0

if save_jci:
   # Save the plot in high res TIFF, new plots will overwrite old ones
   plt.savefig('JCI_plot.tiff', dpi=300, format='tiff')
   print('JCI plot saved as high res TIFF')


# Show the plot
plt.show()
print ()
#in order to organize x y in stored array correctly, need to transpose back again
ssr = np.transpose(ssr)
if f10loop:
    print("With optimized f10 value:", f100min) 
print(f"SSR in F= {round(deltamin, 8)}, r1={round(r1min, 4)} SD={round(sdr1,4)}, r2={round(r2min, 4)} SD={round(sdr2,4)}")
print(f"Theor.SSR={round(rtheor,8)}")
if deltamin>0.3*rtheor:# can use any criterium, can also be done with F-test
    print("Not a good fit, consider re-evaluating the errors in your measurements, zoom in or change weighting scheme")
    
if (r1 > 1 and r2 > 1) or (r1 < 1 and r2 < 1):
    print(f"Azeotrope at f1= {round(1000 * (1 - r2) / (2 - r1 - r2)) / 1000}, "
          "conversion cannot be calculated accurately near that value")
print()


print("f10       weight     f1          Conv     Convcalc     F0         Fexp      Fcalc      delta F    delta Fest")
deltatotal = 0
for sample in range(1, nosample + 1):  # calculate the deltas, etc., for the residuals table
    fnew = f10[sample]
    f2new = 1 - fnew
    F10 = (r1 * fnew ** 2 + fnew * f2new) / (r1 * fnew ** 2 + 2 * fnew * f2new + r2 * f2new ** 2)  # based on instantaneous copolymerization equation
    if Fopt[sample] == f1[sample]:  # avoid divide by zero, etc.
        f1[sample] += 0.00001
    delta = F[sample] - Fopt[sample]
    if abs(delta) > deltaF[sample]:
        remark[sample] += " Error in F too big!, adjust your errors"
    if Fopt[sample] == f1[sample]:
        Fopt[sample] = Fopt[sample] + 0.00001
    deltatotal += w[sample] * (delta) ** 2
    convber = (f10[sample] - f1[sample]) / (Fopt[sample] - f1[sample])
    print(f"{round(f10[sample], 4)}     {round(w[sample], 3)}       {round(f1[sample], 4)}      {round(conv[sample], 3)}      "
          f"{round(convber, 3)}      {round(F10, 3)}      {round(F[sample], 4)}      {round(Fopt[sample], 4)}      "
          f"{round(delta, 4)}      {round(deltaF[sample], 4)}      {remark[sample]}")
print()

print(f"Looked for an optimum of r1 between: {round(r1or - r1range,4)} and {round(r1or + r1range,4)}")
print(f"Looked for an optimum of r2 between: {round(r2or - r2range,4)} and {round(r2or + r2range,4)}")
print()
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import griddata

# Plotting the residuals in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Your existing data and plot setup
start_sample = 1
end_sample = nosample + 1
plot_colors = ['red' if Fopt[i] - F[i] < 0 else 'blue' for i in range(0, nosample + 1)]
ax.scatter(f10[start_sample:end_sample], conv[start_sample:end_sample], F[start_sample:end_sample], c=plot_colors[start_sample:end_sample], marker='o', label='F measured')

# Create a meshgrid from the data ranges
x_range = np.linspace(min(f10[start_sample:end_sample]), max(f10[start_sample:end_sample]), 100)
y_range = np.linspace(min(conv[start_sample:end_sample]), max(conv[start_sample:end_sample]), 100)
X, Y = np.meshgrid(x_range, y_range)

# Interpolate Fopt values on the meshgrid, can be linear or cubic etc
Z = griddata((f10, conv), Fopt, (X, Y), method='cubic')

# Plotting the interpolated plane, does not work very well for strongly curved planes
ax.plot_surface(X, Y, Z, alpha=0.1, color='green')

# Labels and titles
plt.xlabel('$f_{10}$')
ax.set_ylabel('X')
ax.set_zlabel('F')
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.title('3D residuals plot, blue is negative residue, red is positive')

# Saving and displaying

if save_residuals:
   plt.savefig('Residuals_plot.TIFF', dpi=300, format='TIFF')
   print('Residuals plot saved as high-res TIFF')


plt.show()
print()

import csv
 
if save_ssr:
    xyname = input("Give the name of the CSV file to write the sum of squares of residuals space: ")
    while not xyname:
        xyname = input("Please provide a valid name for the CSV file: ")
    xyname = xyname + ".csv"

    with open(xyname, "w", newline='') as file:
        # Specify the delimiter as a comma
        csv_writer = csv.writer(file, delimiter=',')

        # Write header
        csv_writer.writerow([""] + [f"{r1or - r1range + i * stepsizer1}" for i in range(0, 100)])

        # Write data
        for j in range(0, 100):
            r2_value = r2or - r2range + j * stepsizer2
            row_data = [f"{r2_value}"] + [f"{ssr[i, j]}" for i in range(0, 100)]
            csv_writer.writerow(row_data)

    print(f"Data has been written to {xyname}")
    print("r1 on x-axis, r2 on y-axis")



def restart_program():
    root.destroy()  # Destroy the current Tkinter instance
    get_inputs()  # Restart the program


