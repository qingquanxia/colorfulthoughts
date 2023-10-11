---
title: "Unifying the continuous and discrete worlds: the dirac comb"
date: 2023-08-28
draft: false
---

This article is part of a series:
1. [Part 1, The continuous world: Fourier series to Fourier transforms]({{< relref "continuous-fourier.md" >}})
2. [Part 2, Discrete time]({{< relref "discrete-time.md" >}})
3. [Part 3, Unifying the continuous and discrete worlds: the dirac comb]({{< relref "dirac-comb.md" >}})

The Dirac comb is defined as $d_c(t) = \\sum_{n=-\\infty}^{\\infty} \\delta(t-n)$. It is an infinite train of impulses, each spaced exactly one unit apart. To make it spaced $\\Delta t$ apart, we simply use $d_c(t/\\Delta t)$. Notice how this not only adjusts the spacing, but also adjusts the intensity of the impulses, so that on average, integrating over a unit interval will always give you 1. This is a nice property which is analogous to firing paintball pellets at a cart. You can fire slowly with heavy powerful pellets, or you can fire rapidly with small wimpy pellets, but the rate at which momentum is imparted to the cart is the same. In fact, in the limit $\\Delta t\\rightarrow 0$, we have $d_c(t/\\Delta t)\\rightarrow 1$, analogous to hosing the cart down with a continuous stream of water instead of shooting it with discrete pellets.

The Fourier transform of $d_c(t)$ turns out to be itself, that is, $d_c(f)$.

Proof:
Consider first the Fourier series representation of $\\delta(t)$ in the interval $[-1/2,1/2]$. All its Fourier coefficients are equal to 1, hence, in that interval, $\\delta(t) = \\sum_n e^{2\\pi int}$. Because all these exponentials are at least periodic with period 1, we see that $d_c(t) = \\sum_n \\delta(t-n) = \\sum_n e^{2\\pi int}$. Hence, its Fourier trnasform is $FT[d_c(t)] = \\int_t e^{-2\\pi ift} \\sum_n \\delta(t-n) dt = \\sum_n e^{-2\\pi ifn} = d_c(f)$.

The neatest part is that it can be used to model sampling. If we have a continuous signal $x(t)$, and we sample it with spacing $\\Delta t$, what we're really doing is multiplying it by $d_c(t/\\Delta t)$

If we have a signal $x(t)$ contained within a width $T$ and we take the convolution $x(t)*d_c(t/T)$, we are copy-pasting $x(t)$ at every spike in the train, and essentially making a periodic version of the signal with period $T$.

## Deriving the DFT

Suppose we have a sequence $x_n$ which is periodic with period $N$, and spaced $\\Delta t$ apart, hence in the frequency domain it's spaced $\\Delta f = \\frac{1}{N\\Delta t}$ apart. The time domain representation of this sequence is $x(t) = \\sum_n x_n \\delta(\\frac{t}{\\Delta t} - n) = \\sum_n x_n \\Delta t \\delta(t-n\\Delta t)$.

Using linearity, its Fourier transform is
$$\\begin{align}
X(f) &= \\sum_n x_n\\Delta t e^{-2\\pi in\\Delta t f} \\\\
&= \\sum_m \\left( \\sum_{n=0}^{N-1} x_n\\Delta t e^{-2\\pi in\\Delta t f} \\right) e^{-2\\pi imN\\Delta t f} \\\\
&= \\left( \\sum_{n=0}^{N-1} x_n \\frac{1}{N\\Delta f} e^{-2\\pi i \\frac{n}{N}\\frac{f}{\\Delta f}} \\right) \\left( \\sum_m e^{-2\\pi im\\frac{f}{\\Delta f}} \\right) \\\\
&= \\left( \\sum_{n=0}^{N-1} x_n \\frac{1}{N} e^{-2\\pi i \\frac{n}{N}\\frac{f}{\\Delta f}} \\right) \\frac{1}{\\Delta f} d_c(\\frac{f}{\\Delta f}) \\\\
&= \\sum_k X_k \\delta(f-k\\Delta f)
\\end{align}$$

## Applications

#### Sampling

Armed with these facts, proving the Nyquist-Shannon sampling theorem becomes trivial. Suppose you start with a signal $x(t)$. Let's suppose its Fourier transform $X(f)$ is band-limited, that is, contained within an interval of width $F$. Sampling with rate $F$ in the time domain is equivalent to periodizing with period $F$ in the frequency domain. This operation is easily reversed by multiplying by the window function of width $F$ in the frequency domain, which corresponds to convolution by sinc in the time-domain, perfectly reconstructing the original continuous signal.

![](/assets/dirac-comb-sampling.svg)

#### Digital-to-analog converters (DAC)
[TODO: figure sinc overlay of freq]
Suppose we have a digital signal $x(t) = \\sum_{n=-\\infty}^\\infty x[n]\\delta(t/\\Delta t - n)$. When we're ready to output, the DAC uses zero-hold interpolation to reconstruct an analog signal. This is equivalent to convolving the digital signal by the window function $w_{0,\\Delta t}(t)$. Of course, this convolution turns into multiplication by $\\mathrm{sinc}(\\Delta t f) e^{-\\pi i \\Delta t f}$ in the frequency domain. Hence even though it is theoretically possible to perfectly reconstruct a band-limited continuous time signal, in practice the DAC will introduce a distortion. To deal with this distortion, it's important both to choose the operating frequencies of DACs carefully and to place appropriate equalizing filters.


#### Shifting fractional frequency in discrete frequency domain.

Suppose we have a regularly spaced grid of frequency samples, $X_k$, spaced $\\Delta f$ apart, that we have stored in a computer. Shifting by an integer number of spaces is easy. But what if we want to shift by a fraction of a sample, say, $f_s$? The operations we must do is as follows:
1. Convert to time-domain. This results in a signal with period $1/\\Delta f$.
2. Multiply by complex exponential $e^{2\\pi i f_s t}$
3. Multiply by the window function $w_{-\\frac{1}{2\\Delta f},\\frac{1}{2\\Delta f}}(t)$
4. Make it periodic again by convolving by $\\sum_n \\delta(t - \\frac{n}{\\Delta f}) = \\Delta f d_c(\\Delta f t)$
5. Convert back to frequency domain

If we were to do all this in the frequency domain, it would look like:
1. Convolve by $\\delta(f-f_s)$
2. Convolve by $\\frac{1}{\\Delta f}\\mathrm{sinc}(\\frac{f}{\\Delta f})$
3. Multiply by $d_c(\\frac{f}{\\Delta f})$.

Because convolution is associative, we can combine steps 1 and 2 so that it's a single step, that is, convolve by $\\frac{1}{\\Delta f}\\mathrm{sinc}(\\frac{f-f_s}{\\Delta f})$. Step 3 essentially samples on the same grid, and then to get the final coefficients, we need to multiply by $\\Delta f$.

If our original sequence $X_k$ was aperiodic, the shifted frequencies we eventually get is
$$
\\begin{align}
Y_k &= \\left( \\sum_l X_l\\delta(f-l\\Delta f) \\right) * \\mathrm{sinc}(\\frac{f-f_s}{\\Delta f}) \\\\
&= \\sum_l X_l \\mathrm{sinc}(\\frac{f-f_s-(l\\Delta f)}{\\Delta f})\\\\
&= \\sum_l X_l \\mathrm{sinc}(k-l-\\frac{f_s}{\\Delta f})
\\end{align}
$$




#### TODO: Concatenating OFDM symbols, and windowing
In communication systems, it is common to divide the infinite line of time into so called "symbols". For example, in 4G LTE, a symbol is about 70us. Due to the nature of the Fourier series, each of these symbols can be represented by countably many "subcarriers" regularly spaced in a frequency lattice. When these symbols are concatenated together, they are often discontinuous at the boundaries, which results in out-of-band emissions. The reason is intuitively clear. In the time domain, doing such a thing is essentially taking your periodic signal, which represents one symbol, and multiplying it by the window function. In the frequency domain, this is convoluting by a sinc function. The smaller your symbol is, the more broad the sinc function, and the more out-of-band emissions you'll produce. However, if you smooth the window out at its boundaries, the function you convolute by in the frequency domain is more localized.


## TODO: Power

The energy of a signal in some interval $[a,b]$ is defined as $\\int_{a}^{b} x^* (t) x(t) dt$. The total energy, $\\int_{t=-\\infty}^{\\infty} x^* (t) x(t) = \\braket{x|x} = |x|^2$, is a property of the signal itself, and therefore stays the same if we change basis into frequency domain. Hence, $\\int_{t=-\\infty}^{\\infty} x^* (t) x(t) = \\int_{f=-\\infty}^{\\infty} X^* (f) X(f)$. This is called Parseval's theorem.

Unfortunately, energy doesn't play so well with sampled signals, because squaring a dirac delta yields something which integrates to infinity. One trick we can play, though, is if we have a signal $x(t)$, and we sample it to get $x(t)d_c(t/\\Delta t)$, then an approximation of the original energy is $\\braket{x|x} \\approx \\braket{x(t)d_c(t/\\Delta t) | x(t)}$.

The instantaneous power at time $t$ is defined as the rate of energy delivered, that is, $x^* (t) x(t)$. With this definition, any pure tone $e^{2\\pi ift}$ has power 1.

Suppose we wanted to plot some sort of power spectrum on a spectrum analyzer. Our instrument is finite, so let's say it captures the signal in the interval $[0,T]$, and then extends the signal to infinity by periodizing. Then it gets a frequency domain plot with spacing $\\Delta f = 1/T$. The total power is simply $\\sum_n |X(n)|^2$. If we doubled the interval, the frequency spectrum wouldn't change: $\\Delta f$ would decrease by two, but we would now have $X(nf_0)=0$ for all odd $n$. Hence it makes sense to speak of the units as dBm/Hz.
