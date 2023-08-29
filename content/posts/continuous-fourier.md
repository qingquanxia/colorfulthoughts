---
title: "Continuous Fourier"
date: 2023-08-26T16:26:48-07:00
draft: false
---

This series will consist of four parts:
1. Fourier series to the Fourier transform
2. The DFT and DT fourier series
3. Unifying the continuous and discrete worlds: the dirac comb



There is something special about waves. You can see them as the oscillating ocean surface, bobbing ships up and down as they travel to far off places. You can see them when you strike a long metal hand-rail, as the metal vibrates and causes the pole to be blurry in your vision, and then hear them as the metal transfers its vibration into the pressure of the air. You can see them quite literally, because light is just a ripple in the invisible electromagnetic field vibrating at $10^{15}$ cycles per second. The purpose of this page is to explain what a wave is in mathematical terms.

The exponential function:

$e^z$ is defined as a power series $1+z+\\frac{z^2}{2!}+\\frac{z^3}{3!}+\\cdots$. It is a fact that if you plot $e^{i\\theta}$ in the complex plane, as you increase $\\theta$ you will see a particle starting from the x-axis, rotating counterclockwise in a circle at constant speed. From this we define sin and cos as the real and imaginary parts: $\\cos(\\theta) := \\Re{e^{i\\theta}} = (e^{i\\theta}+e^{-i\\theta})/2$, and $\\sin(\\theta) := \\Im{e^{i\\theta}} = (e^{i\\theta}-e^{-i\\theta})/(2i)$. Euler's formula summarizes this decomposition nicely: $e^{i\\theta} = \\cos(\\theta) + i\\sin(\\theta)$.

Differentiating the power series reveals something interesting: $d/dx(e^x) = e^x$ (I am restricting myself to functions $\\RR\\rightarrow\\CC$. This also holds true for $\\CC\\rightarrow\\CC$, but that requires complex analysis to explain, and isn't really relevant right now). This gives us a really neat alternative definition: "The exponential function is defined as the unique eigenfunction (up to scale factor) of the differentiation operator".

Because it turns out pretty much every formula in physics is a differential equation, the operator $(d/dt)$ is extremely important. That is why the exponential function has been called the single most important function in math. In general the derivative does messy things to your function. But if you can express your function as the sum of exponentials, now you've got some serious simplification power.

Inner product spaces and Dirac notation

The idea of inner products originates from the 3D Euclidean vectors used in classical mechanics. A vector $v$ is an arrow in 3D space; if you place the tail of the arrow at the origin, you can associate the vector with the coordinates at the head, $(v_1, v_2, v_3)$. The dot product of two vectors $u,v$ is defined as $u\\cdot v = u_1v_1+u_2v_2+u_3v_3$. It essentially measures how aligned they are. Two vectors that are perfectly aligned will give maximum dot product; perpendicular will give 0; anti-parallel will give most negative dot product. Notice in particular that $v\\cdot v = |v|^2$; in fact this property can be used to define the length of a vector. They are also used for decomposing a vector into coordinates. If I have an orthonormal basis $\\{e_1, e_2, e_3\\}$, I can project my vector $v$ onto the lines spanned by each of them, so that projecting it onto $e_1$ will give the vector $(v\\cdot e_1)e_1$. Because this orthonormal basis spans the entire 3D space, adding the three projection vectors will perfectly reconstruct the original vector.

Of course, people just weren't content with sticking to three dimensions. The same concepts can be applied to arbitrary dimension, so that for $u,v\\in \\RR^n$, $u\\cdot v = \\sum_i u_iv_i$. The real fun begins when we try to extend this concept to complex vectors, that is, to $\\CC^n$. Because length is absolutely essential to geometry, we would like to preserve the idea that the product of a vector with itself is its length squared, that is, it is real. Hence we will define the product of two vector $u,v\\in \\CC^n$ as $\\sum_i u_i^* v_i$. Hence, the product of $v$ with itself is $\\sum_i v_i^* v_i = \\sum_i |v_i|^2$, which makes sense.

Here we formalize the concept of an inner product, which is a generalization of all the ideas expressed above. An inner product over a complex (or real) vector space $V$ is a product $V\\times V\\R \\CC$ which has these three properties:
1. Linearity: $\\braket{u|a v + b w} = a\\braket{u|v} + b\\braket{u|w}$
2. Conjugate symmetry: $\\braket{u|v} = \\braket{v|u}^*$
3. Positive definite: $\\braket{v|v}$ is greater than zero if $v\\neq 0$, and $0$ otherwise. In particular, its value is real.

A few remarks are in order. First of all, this definition works equally well with a real vector space, since the complex conjugate of a real number is just itself. Second, properties 1 and 2 combine to give $\\braket{au + bv| w} = a^* \\braket{u|w} + b^* \\braket{v|w}$, which can be tricky to remember. Third, these definitions seem more general than our previous complex-coefficient based way of calculating inner products, but it is actually not really. If we are given a finite dimensional complex inner product space $V$, we are guaranteed by the Gram-Schmidt procedure to find an orthonormal basis $\\{e_1,\\cdots,e_n\\}$, where $\\braket{e_i|e_j}$ is 0 if $i\\neq j$ and 1 if $i=j$. Hence the inner product of any two vectors is, by linearity, $\\braket{u|v} = \\braket{\\sum_i u_i e_i | \\sum_j v_j e_j} = \\sum_{i,j} u_i^* v_j \\braket{e_i | e_j} = \\sum_i u_i^* v_i$, which is precisely our old definition.

(An aside on notation... Quantum mechanics has to deal with this kind of math so much that Paul Dirac, one of the early scientists in the field, decided to take the notation further. He called vectors "kets", and wrote them like this: $\\ket{v}$. Each vector has a dual-vector called a "bra", which is the conjugate transpose: $\\bra{v}$. Hence the matrix product of a bra and ket corresponded exactly to the inner product of the two vectors, which he called a "bra-ket": $\\braket{u|v}$. Very funny Dirac. You can combine them in various ways. For example, if $A,B$ are linear operators, then $\\braket{Au|B|v} = u^\\dag A^\\dag B v$.)

So our new definition is still the same familiar product. However, it is independent of basis, and it does give us a new way to apply inner products to the strange world of uncountably-infinite dimensional vector spaces...

Fourier stuff.

To keep the math consistent, in the pure mathy parts I will always use $t$ as the independent variable, $x,y$ as functions of $t$, and $f$ to denote frequency. Keep in mind the same math can be applied whether your independent variable is space or time. In the physics parts however I reserve the right to go crazy with variable names.

Consider the space of continuous functions $[a,b]\\rightarrow \\CC$. It is a vector space, and we can define an inner product on it: $\\braket{x(t)|y(t)} := \\int_a^b x^*(t)y(t)dt$. Now let's consider what the inner product of two exponentials looks like: $\\braket{e^{2\\pi if_1 t}| e^{2\\pi if_2 t}} = \\int_a^b e^{2\\pi i(f_2-f_1)t}dt = [\\frac{1}{2\\pi i(f_2-f_1)}e^{2\\pi i(f_2-f_1)}]_a^b$. In order for this to be zero, we want $e^{2\\pi i(f_2-f_1)b} = e^{2\\pi i(f_2-f_1)a}$, that is, $2\\pi (f_2-f_1)(b-a) = 2\\pi n$ for some integer $n$. So a necessary and sufficient condition for two exponentials to be orthogonal is that $f_2-f_1 = n/(b-a)$ for some integer $n$.

This leads us to consider the set of functions $\\{ \\frac{1}{\\sqrt{b-a}} e^{2\\pi ifx} \\}_{f=n/(b-a)}$ for all integers $n$. It's definitely an orthonormal set of functions. Is it complete? You'll have to take my word for it, but it is in fact complete. (Or, if the idea of trusting a complete stranger from who knows where with a result as important as this bothers you, you should take a look at a book like "Mathematics of Classical and Quantum Physics" by Byron and Fuller, which is filled with wonderful proofs of these kinds of things).

We have now decomposed any function on the interval into the sum of complex exponentials! This lets us do interesting things....

Physics aside: Guitar strings
Consider a string stretched tight between two fixed points. Suppose its mass per unit length is $\\lambda$, and it is subject to a tension of $T$.


From the diagram above, we see that the y-component of force on the tiny bit of string is $-T\\frac{dy}{dx} + T(\\frac{dy}{dx} + \\frac{d^2y}{dx^2}dx) = T\\frac{d^2y}{dx^2}dx$. Hence, by Newton's second law $F=ma$, we have $T\\frac{d^2y}{dx^2}dx = (\\lambda dx) (\\frac{d^2y}{dt^2})$, which gives us a "wave equation": $\\frac{d^2y}{dx^2} = \\frac{1}{T/\\lambda} \\frac{d^2y}{dt^2}$.

For comparison, here are some other wave equations:
1. For sound waves: $\\nabla^2 p = \\frac{1}{v_s^2}\\frac{\\partial^2 \\bf{E}}{\\partial t^2}$ where $v_s=\\sqrt{B/\\rho_0}$ is the speed of sound. $B$ is the bulk modulus of the fluid, and $\\rho_0$ is the average density.
2. For EM waves (light): $\\nabla^2 \\bf{E} = \\frac{1}{c^2}\\frac{\\partial^2 \\bf{E}}{\\partial t^2}$, where $c=1/\\sqrt{\\epsilon_0\\mu_0}$ is the speed of light.
3. For matter waves (like free electrons) in quantum theory: $\\nabla^2 \\psi = \\frac{2m}{\\hbar} (-i\\frac{\\partial \\psi}{\\partial t})$. (This looks a bit different from the others; it expresses the fact that matter waves do not have a constant phase-velocity, which leads to dispersion of different frequency components. But keep in mind the $i$ is essentially a derivative, so it has sorta the same form).
4. For gravitational waves: I have no idea, but in the limit of small distortions I'd guess it has the same form as the others.

(These equations are in contrast with the heat and diffusion equations, in which we only differentiate with respect to time once. Hence these do not exhibit wave-like behavior.)

With foresight, I will define $v=\\sqrt{T/\\lambda}$.

To solve this equation, let's start by only looking for solutions of the form $y(x,t) = f(x) g(t)$. This seems a ridiculously strict assumption at first, but we will see soon that it actually yields very general solutions.

Plugging it in, we get $\\delsq{f}{x} g = \\frac{1}{v^2}f\\delsq{g}{t}$ so $\\frac{1}{f}\\delsq{f}{x} = \\frac{1}{v^2}\\frac{1}{g}\\delsq{g}{t}$. Now here is a subtle point. The right side of the equation depends only on $t$, while the right side of the equation depends only on $x$. This means that both sides are equal to some constant, which, having foresight, I will name $-k^2$. If it were not so, then I could change the value of the left hand side by moving $t$ and keeping $x$ constant; because of the equality, the right side will also be changed; but this is a contradiction, because the right side depended only on $x$.

Now we have essentially decomposed our PDE into two ODEs, which are simple to solve. We end up with the general solutions $f(x) = Ae^{ikx} + Be^{-ikx}$ and $g(t) = Ce^{i\\omega t} + De^{-i\\omega t}$. The boundaries conditions specify that both ends of the string are fixed at zero, hence $A+B=0$ and $Ae^{ikL}+Be^{-ikL}=0$. Substituting the first equation into the second gives $A(e^{ikL}-e^{-ikL})=0 \\Rightarrow e^{ikL}=e^{-ikL}\\Rightarrow e^{2ikL}=1$, which implies $k=n\\pi/L$ for some integer $n$. Taking into account that the wave is real, we get as our final solution: $y(x,t) = \\sin(kx)e^{i\\omega t} + \\text{c.c.}$. The c.c. stands for complex-conjugate, and it is just there to make the solution real. Some people will set in place a convention that the final solution is complex, and it's implicit that the physical solution is obtained by taking the real part. This is totally fine for linear systems, and is commonly used in discussing phasors in electrical circuits, however it gets a bit confusing when non-linearity comes into play, because $\\Re{ab} \\neq \\Re{a}\\Re{b}$.

The solutions we have found are called the standing waves, because they oscillate in place. They are the first, second, third, ... harmonics, all multiples of the fundamental harmonic. They seem like special cases, until you realize that using Fourier series, we can decompose ANY function into a sum of these standing waves.


Now let's say we pluck the string at $x=0.35L$, so it looks like a triangle at $t=0$.





Fourier series can be used to model all sorts of interesting boundary problems, from the electromagnetic field in your microwave (metal walls impose zero electric field boundary conditions), to the variation of sound in an accoustic concert hall, to the behavior of an electron trapped in an infinite potential well.

They run into a limitation, however. Time is essentially infinite in both directions. Space is also infinite in all directions. How do we model light waves from distant stars, which took years to reach us?

We need to consider the limit as $[a,b]$ stretches out into an infinite interval.

Let $a\\rightarrow -\\infty$ and $b\\rightarrow +\\infty$. Our inner product becomes $\\braket{x(t),y(t)} \\approx \\int_{-\\infty}^{+\\infty}x^*(t)y(t)dt$. We can still write an arbitrary function $x(t)$ as $x(t) = \\sum_{n=-\\infty}^\\infty  \\braket{\\frac{1}{\\sqrt{b-a}}e^{2\\pi i\\frac{n}{b-a}t'}, x(t')} \\frac{1}{\\sqrt{b-a}}e^{2\\pi i\\frac{n}{b-a}t} \\approx \\int_{f=-\\infty}^\\infty \\braket{e^{2\\pi ift'}, x(t')} e^{2\\pi ift} df = \\int_{f=-\\infty}^\\infty [\\int_{t'=-\\infty}^\\infty e^{-2\\pi ift'} x(t') dt'] e^{2\\pi ift} df$.

There we have it! We have derived our Fourier transform:
$X(f) = \\int_{-\\infty}^\\infty e^{-2\\pi ift} x(t)dt$
$x(t) = \\int_{-\\infty}^\\infty e^{2\\pi ift} X(f)df$

These equations are remarkably symmetric. In fact, it turns out that $FT[x(t)] = -IFT[x(-t)]$.

This has its uses. For example, one question that the Fourier transforms immediately raise is: what is the Fourier transform of $e^{2\\pi if_0t}$? If we try plugging it in, it yields a strange integral. For any $f\\neq f_0$, that part of the integral should be zero, since a complex exponential is, over the whole real line, zero on average. However, when $f=f_0$, the integral is infinity. This kinda makes sense since a pure tone should yield an infinitely sharp response in the frequency domain.

But we can be more precise. What is the Fourier transform of a dirac delta? $FT[\\delta (t-t_0)] = \\int_{t=-\\infty}^\\infty e^{-2\\pi ift} \\delta(t-t_0)dt = e^{-2\\pi ift_0}$. But from duality, this means $FT[e^{2\\pi if_0t}] = \\delta(f-f_0)$!

A different perspective: quantum mechanics and momentum space
I first learned about the Fourier transform in quantum mechanics, but it was presented in such a way that I did not even know it was a transform. Its description was so natural and intuitive that it just made sense. What follows is a quick debrief on quantum mechanics, but the math can be applied to suprisingly many other physical systems.
There are five postulates in quantum mechanics:
1. A physical system's state is a vector in a complex Hilbert space. (A Hilbert space is just an inner-product space whose norm makes it a complete metric space. All the familiar finite dimensional spaces are Hilbert spaces. $L^2(\\RR, \\CC)$ is also a Hilbert space).
2a. A physical observable (e.g. energy, momentum, position, spin) is represented by a Hermitian linear operator on the space. (An operator $A$ is Hermitian if its adjoint, aka its conjugate-transpose, is equal to itself.) When you make a measurement of the observable using any sort of apparatus, it is only possible to observe an eigenvalue of the Hermitian operator.
2b. [not relevant to us] Measurement probability postulate.
2c. [not relevant to us] Collapse postulate.
3. [not relevant to us] Time evolution.

For a particle that is free to travel along the entire 1-dimensional line, the Hilbert space we are interested in is $L^2(\\RR, \\CC)$.

There are two operators we are interested in. The position observable is given by the operator $f(x) \\mapsto xf(x)$, that is, multiplication by $x$. It is easy to check that this is linear and Hermitian. Its eigenvalues consist of the entire real line. An eigenvalue $c$ has eigenvector $\\delta(x-c)$. Because its eigenvectors form an orthonormal basis for the Hilbert space, we can express any wave $\\psi(x)$ in its basis. In this case it is trivial: $\\ket{\\psi} = \\int_{-\\infty}^{+\\infty} \\psi(x) \\ket{\\delta(x-c)} dx$.

The momentum observable, on the other hand, is the operator $\\hbar/i (d/dx)$. Again, this is both linear and Hermitian. Its eigenvectors comprise of the solution of the differential equation $\\hbar/i (d/dx)\\phi(x) = p\\phi(x)$, which turn out to be the exponential functions!

So we have two different sets of orthonormal bases for our Hilbert space. Getting from one from the other is a simple matter of changing bases. The operator $\\int_p \\ket{p}\\bra{p} dp$ is called a "resolution of the identity", because it is in fact the identity; it merely breaks our vector up as a linear combination of the basis. If we expand it out, what it is really saying is that $\\ket{\\psi} = (\\int_p \\ket{p}\\bra{p} dp)\\ket{\\psi} = \\int_p (^{ipx/\\hbar} (\\braket{p|\\psi})dp = \\int_p e^{ipx/\\hbar} (\\int_{x'} e^{-ipx'/\\hbar}\\psi dx')dp$. But this is just our FT/IFT! So in essence, the Fourier transform is nothing more than a change in basis.










The DFT and Discrete Fourier Series.
Because computers can only deal with finite quantities, any signal processing done inside a computer must be finite.
Following the example of Fourier series, we now define the discrete fourier transform:

$DFT[x(n)] =\\sum{n=0}^{N-1} e^{-2\\pi ifn/N} x(n)$
$IDFT[x(n)] =\\frac{1}{N}\\sum{n=0}^{N-1} e^{-2\\pi in/N} x(n)$





Unifying the continuous and discrete worlds: the dirac comb
The Dirac comb is defined as $d_c(t) = \\sum_{n=-\\infty}^{\\infty} \\delta(t-n)$. It is an infinite train of impulses, spaced exactly one unit apart. To make it spaced $T$ apart, we simply use $d_c(t/T)$. Notice how this not only adjusts the spacing, but also adjusts the intensity of the impulses so that on average, integrating over a unit interval will give you 1. This is a nice property which is analogous to firing paintball pellets at a cart. You can fire slowly with heavy powerful pellets, or you can fire rapidly with small wimpy pellets, but the rate at which momentum is imparted to the cart is the same. In fact, in the limit $T\\rightarrow 0$, we have $d_c(t/T)\\rightarrow 1$, analogous to hosing the cart down with a continuous stream of water instead of shooting it with discrete pellets.

The Fourier transform of $d_c(t)$ is actually itself, that is, $d_c(f)$.

The neatest thing about it is it can be used to model sampling. If we have a continuous signal $x(t)$, and we sample it with period $T$, which we're really doing is multiplying it by $d_c(t/T)$

If we have a signal $x(t)$ contained within a width $T$ and we take the convolution $x(t)*d_c(t/T)$, we are copy-pasting $x(t)$ at every spike in the train, and essentially making a periodic version of the signal with period $T$. (If the idea of convoluting by an infinite train of stuff bothers you, keep in mind that if our signal is square-integrable, then its square typically decreases faster than $1/x^2$, so the sum will converge at every point since $\\sum_{n=0}^\\infty 1/n^2$ is finite.)

[figure Nyquist]

With these facts, proving the Nyquist sampling theorem becomes trivial. Suppose you start with a signal $x(t)$. Let's suppose its Fourier transform $X(f)$ is band-limited, that is, contained within an interval of width $f_b$. Sampling with rate $f_b$ in the time domain is equivalent to periodizing with period $f_b$ in the frequency domain. This operation is easily reversed by multiplying by the window function of width $f_c$ in the frequency domain, which corresponds to convolution by sinc in the time-domain, perfectly reconstructing the original continuous signal.

Digital-to-analog converters (DAC)
[figure sinc overlay of freq]
Suppose we have a continuous valued signal $x(t)$ which we feed into an ADC. This ADC samples it with period $T$, so in the digital domain we have the signal $x(t)d_c(t/T)$. When it comes out, the DAC uses zero-hold interpolation to reconstruct an analog signal. This is equivalent to convoluting the digital signal by the window function $w_{0,T}(t)$. Of course, this convolution in the time domain turns into multiplication by $\\sinc(fT)$ in the frequency domain. Hence even though it is theoretically possible to perfectly reconstruct a band-limited continuous time signal, in practice the DAC will introduce a distortion. This is why it is very important to choose the operating frequencies of DACs carefully, to avoid significant distortion.


Shifting fractional frequency in discrete frequency domain.

Suppose we have a regularly spaced grid of frequency samples, spaced $f_0$ apart, that we have stored in a computer. Shifting by an integer number of spaces is easy. But what if we want to shift by a fraction of a sample, say, $\\Delta f$? The operations we must do is as follows:
1. Convert to time-domain. This results in a signal with period $1/f_0$.
2. Multiply by complex exponential $e^{2\\pi i \\Delta f t}$
3. Multiply by the window function $w_{0,1/f_0}(t)$
4. Make it periodic again by convoluting by $d_c(f_0 t)$
5. Convert back to frequency domain

If we were to do all this in the frequency domain, it would look like:
2. Convolve by $\\delta(f-\\Delta f)$
3. Convolve by $\\sinc(...)$
4. Sample again in the same grid.

The formula we eventually get is something like this:
$X'(f_0 n) = \\sum_m \\sinc(...) X'(f_0 m)$




Concatenating OFDM symbols, and windowing
In communication systems, it is common to divide the infinite line of time into so called "symbols". Due to the nature of the Fourier series, each of these symbols has countably many "subcarriers" regularly spaced in a frequency lattice. When these symbols are concatenated together, they are often discontinuous at the boundaries, which results in out-of-band emissions. The reason is intuitively clear. In the time domain, doing such a thing is essentially taking your periodic signal, which represents one symbol, and multiplying it by the window function. In the frequency domain, this is convoluting by a sinc function. The smaller your symbol is, the more broad the sinc function, and the more out-of-band emissions you'll produce. However, if you smooth the window out at its boundaries, the function you convolute by in the frequency domain is more localized.


Power
The energy of a signal in some interval $[a,b]$ is defined as $\\int_{a}^{b} x(t)^* x(t) dt$. The total energy, $\\int_{t=-\\infty}^{\\infty} x(t)^* x(t) = \\braket{x,x} = |x|^2$, is a property of the signal itself, and therefore stays the same if we change basis into frequency domain. Hence, $\\int_{t=-\\infty}^{\\infty} x(t)^* x(t) = \\int_{f=-\\infty}^{\\infty} X(f)^* X(f)$. This is called Parseval's theorem.

Unfortunately, energy doesn't play so well with sampled signals, because squaring a dirac delta yields something which integrates to infinity. One trick we can play, though, is if we have a signal $x(t)$, and we sample it to get $x(t)d_c(t/T)$, then an approximation of the original energy is $\\braket{x|x} \\approx \\braket{x(t)d_c(t/T) | x(t)}$.

If the signal is long and approximately periodic, it makes sense to define the average power as the average amount of energy delivered per unit time, that is, $\\frac{1}{T} \\int_a^{a+T} x(t)^* x(t) dt$. With this definition, any pure tone $e^{2\\pi ift}$ has power 1.

Suppose we wanted to plot some sort of power spectrum on a spectrum analyzer. Our instrument is finite, so let's say it captures the signal in the interval $[0,T]$, and then extends the signal to infinity by periodizing. Then it gets a frequency domain plot with spacing $f_0 = 1/T$. The total power is simply $\\sum_n |X(n)|^2$. If we doubled the interval, the frequency spectrum wouldn't change: $f_0$ would decrease by two, but we would now have $X(nf_0)=0$ for all odd $n$. Hence it makes sense to speak of the units as dBm/Hz.
