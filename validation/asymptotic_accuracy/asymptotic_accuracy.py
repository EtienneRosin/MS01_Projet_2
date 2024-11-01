import numpy as np
import matplotlib.pyplot as plt

style_path = "validation/utils/custom_science.mplstyle"

jacobi_file_name = "validation/asymptotic_accuracy/jacobi_asymptotic_acuracy_measurement.csv"
conjugate_gradient_file_name = "validation/asymptotic_accuracy/conjugate_gradient_asymptotic_acuracy_measurement.csv"
file_props = dict(dtype = float, delimiter = ',', names=True)

jacobi_data = np.genfromtxt(fname = jacobi_file_name, **file_props)
conjugate_gradient_data = np.genfromtxt(fname = conjugate_gradient_file_name, **file_props)
# print(data['h'])
with plt.style.context(style_path):
    fig = plt.figure()
    ax = fig.add_subplot()

    line_props = dict(marker = "o", linestyle = "--", linewidth = 1)

    ax.plot(np.log(jacobi_data['h']), np.log(jacobi_data['error']), label = "Jacobi", **line_props)
    ax.plot(np.log(conjugate_gradient_data['h']), np.log(conjugate_gradient_data['error']), label = "C. G.", **line_props, markersize = 3)
    ax.plot(np.log(jacobi_data['h']), np.log(jacobi_data['h']**2), label = r"$h \mapsto h^2$", linestyle = "--")
    ax.set(xlabel = r"$\log(h)$", ylabel = r"$\log(\|e_h\|)$")
    ax.legend()
    fig.savefig(fname="validation/asymptotic_accuracy/asymptotic_accuracy.pdf")
    fig.savefig(fname="validation/asymptotic_accuracy/asymptotic_accuracy.svg")
    plt.show()