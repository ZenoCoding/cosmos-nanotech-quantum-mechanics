# Animates a function f(x) over time t.
# (designed with the time evolution of the diffusion equation in mind,
#  but works for a general function).
#
# Run it like this:
#
#     python3 animatedif.py input.txt xmin xmax tmax
#
# where input.txt (the file extension can be anything) is the text file containing
# the data to animate organized as follows:
#
#     fx00 fx01 fx02 ...
#     fx10 fx11 fx12 ...
#     fx20 fx21 fx22 ...
#     .
#     .
#
# that is, each row of the text file contains the value of f at each spatial
# point (in increasing x) and increasing time.
#
# Also input are xmin and xmax: the smallest x value f is evaluated at (xmin)
# and the largest one (xmax).
# And tmax: the time (from 0) the diffusion equation was simulated from (row 0
# of the text file is f at t = 0, row 1 t = dt, row 2 t = 2dt, ... where dt = tmax / nrows).
#
import sys
import numpy
import matplotlib.pyplot
import matplotlib.animation

# adjust me to make the animation faster/slower
# lower  = faster (min 1)
# higher = slower
interval = 10


def initial(curve, timer):
    curve.set_data([], [])
    timer.set_text("")
    return curve, timer

def f(frame, ts, xs, fxs, curve, timer):
    curve.set_data(xs, fxs[frame, :])
    timer.set_text(f"t = {ts[frame]:.3f}")
    return curve, timer

fname = sys.argv[1]
xmin = float(sys.argv[2])
xmax = float(sys.argv[3])
tmax = float(sys.argv[4])

fxs = numpy.loadtxt(fname)
nts, nxs = fxs.shape

figure, axes = matplotlib.pyplot.subplots(figsize=(8,6), dpi=200)

xdiff = xmax - xmin
ymin = 0.0 ; ymax = numpy.max(fxs)
ydiff = ymax - ymin
axes.set_xlim(xmin - 0.05*xdiff, xmax + 0.05*xdiff)
axes.set_ylim(ymin - 0.05*ydiff, ymax + 0.05*ydiff)

xs = numpy.linspace(xmin, xmax, nxs)
ts = numpy.linspace(0, tmax, nts)


axes.set_xlabel("x", fontsize=10)
axes.set_ylabel("f(x)", fontsize=10)

axes.tick_params(axis="both", which="major", labelsize=10, length=6, width=1)

axes.grid(True)

curve, = axes.plot([], [], lw=2, color="red")
timer = axes.text(x=0.95, y=0.95, s="", transform=axes.transAxes, fontsize=10, ha="right")

animation = matplotlib.animation.FuncAnimation(figure,
                                               f,
                                               frames=nts,
                                               init_func=lambda: initial(curve, timer),
                                               fargs=(ts, xs, fxs, curve, timer),
                                               interval=interval,
                                               blit=True)

# Show the animation
matplotlib.pyplot.show()
# Alternatively you can save it as a file (may take a while without a good "writer"):
# animation.save("animation.gif")