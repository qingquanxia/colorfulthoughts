---
title: "The continuous world: Fourier series to Fourier transforms"
date: 2023-08-26
draft: false
---

This article is part of a series:
1. [Part 1, The continuous world: Fourier series to Fourier transforms]({{< relref "continuous-fourier.md" >}})
2. [Part 2, Discrete time]({{< relref "discrete-time.md" >}})
3. [Part 3, Unifying the continuous and discrete worlds: the dirac comb]({{< relref "dirac-comb.md" >}})

There is something special about waves. You can see them as the oscillating ocean surface, bobbing ships up and down as they travel to far off places. You can see them when you strike a long metal hand-rail, as the metal vibrates and causes the pole to be blurry in your vision, and then hear them as the metal transfers its vibration into the pressure of the air. You can see them quite literally, because light is just a ripple in the invisible electromagnetic field vibrating at $10^{15}$ cycles per second. The purpose of this page is to explain what a wave is in mathematical terms.

## The exponential function

$e^z$ is defined as a power series $1+z+\\frac{z^2}{2!}+\\frac{z^3}{3!}+\\cdots$. It is a fact that if you plot $e^{i\\theta}$ in the complex plane, as you increase $\\theta$ you will see a particle starting from the x-axis, rotating counterclockwise in a unit circle at constant speed. From this we define sin and cos as the real and imaginary parts: $\\cos(\\theta) := \\Re{e^{i\\theta}} = (e^{i\\theta}+e^{-i\\theta})/2$, and $\\sin(\\theta) := \\Im{e^{i\\theta}} = (e^{i\\theta}-e^{-i\\theta})/(2i)$. Euler's formula summarizes this decomposition nicely: $e^{i\\theta} = \\cos(\\theta) + i\\sin(\\theta)$.

Differentiating the power series reveals something interesting: $(d/dt)(e^t) = e^t$ (I am restricting myself to functions $\\RR\\R\\CC$. This also holds true for $\\CC\\R\\CC$, but that requires complex analysis to explain, and isn't really relevant right now). This gives us a really neat alternative definition: "The exponential function is defined as the unique eigenfunction of the differentiation operator with eigenvalue 1, such that $e^0=1$".

Because it turns out pretty much every formula in physics is a differential equation, the operator $(d/dt)$ is extremely important. That is why the exponential function has been dubbed "the most important function in mathematics" (Walter Rudin). In general the derivative does messy things to your function. But if you can express your function as the sum of exponentials, now you've got some serious simplification power.

## Inner product spaces and Dirac notation

The idea of inner products originates from the 3D Euclidean vectors used in classical mechanics. A vector $v$ is an arrow in 3D space; if you place the tail of the arrow at the origin, you can associate the vector with the coordinates at the head, $(v_1, v_2, v_3)$. The dot product of two vectors $u,v$ is defined as $u\\cdot v = u_1v_1+u_2v_2+u_3v_3$. It essentially measures how aligned they are. Two vectors that are perfectly aligned will give maximum dot product; perpendicular will give 0; anti-parallel will give most negative dot product. Notice in particular that $v\\cdot v = |v|^2$; in fact this property can be used to define the length of a vector. They are also used for decomposing a vector into coordinates. If I have an orthonormal basis $\\{e_1, e_2, e_3\\}$, I can project my vector $v$ onto the lines spanned by each of them, so that projecting it onto $e_1$, for example, will give the vector $(v\\cdot e_1)e_1$. Because these three basis vectors span the entire 3-dimensional space, adding the three projection vectors will perfectly reconstruct the original vector.

Of course, people just weren't content with sticking to three dimensions. The same concepts can be applied to arbitrary dimension, so that for $u,v\\in \\RR^n$, $u\\cdot v = \\sum_i u_iv_i$. The real fun begins when we try to extend this concept to complex vectors, that is, to $\\CC^n$. Because length is absolutely essential to geometry, we would like to preserve the idea that the product of a vector with itself is its length squared, that is, it is real. Hence we will define the product of two vector $u,v\\in \\CC^n$ as $\\sum_i u_i^* v_i$. Hence, the product of $v$ with itself is $\\sum_i v_i^* v_i = \\sum_i |v_i|^2$, which makes sense.

Here we formalize the concept of an inner product, which is a generalization of all the ideas expressed above. An inner product over a complex (or real) vector space $V$ is defined as a function mapping $V\\times V\\R \\CC$ which has these three properties:
1. Linearity: $\\braket{u|a v + b w} = a\\braket{u|v} + b\\braket{u|w}$
2. Conjugate symmetry: $\\braket{u|v} = \\braket{v|u}^*$
3. Positive definite: $\\braket{v|v}$ is greater than zero if $v\\neq 0$, and $0$ otherwise. In particular, its value is real.

A few remarks are in order.

First of all, this definition works equally well with a real vector space, since the complex conjugate of a real number is just itself.

Second, properties 1 and 2 combine to give $\\braket{au + bv| w} = a^* \\braket{u|w} + b^* \\braket{v|w}$, which can be tricky to remember.

Third, these definitions seem more general than our previous complex-coefficient based way of calculating inner products, but it is actually not. If we are given a finite dimensional complex inner product space $V$, we are guaranteed by the Gram-Schmidt procedure to find an orthonormal basis $\\{e_1,\\cdots,e_n\\}$, where $\\braket{e_i|e_j}$ is 0 if $i\\neq j$ and 1 if $i=j$. Hence the inner product of any two vectors is, by linearity, $\\braket{u|v} = \\braket{\\sum_i u_i e_i | \\sum_j v_j e_j} = \\sum_{i,j} u_i^* v_j \\braket{e_i | e_j} = \\sum_i u_i^* v_i$, which is precisely our old definition.

(An aside on notation... Quantum mechanics has to deal with this kind of math so much that Paul Dirac, one of the early scientists in the field, decided to take the notation further. He called vectors "kets", and wrote them like this: $\\ket{v}$. Each vector has a dual-vector called a "bra", which is the conjugate transpose: $\\bra{v}$. Hence the matrix product of a bra and ket corresponded exactly to the inner product of the two vectors, which he called a "bra-ket": $\\braket{u|v}$. Very funny Dirac. You can combine them in various ways. For example, if $A,B$ are linear operators, then $\\braket{Au|B|v} = u^\\dagger A^\\dagger B v$.)

So our new definition is still the same familiar product. However, it is independent of basis, and it does give us a new way to apply inner products to the strange world of uncountably-infinite dimensional vector spaces.

## TODO: spectral theorem? Should I put it here or after fourier transforms?

## Fourier series

To keep the math consistent, in the pure mathy parts I will always use $t$ as the independent variable, $x,y$ as functions of $t$, and $f$ to denote frequency. Keep in mind the same math can be applied whether your independent variable is space or time. In the physics parts however I reserve the right to go crazy with variable names.

Consider the space of square-integrable functions $[a,b]\\rightarrow \\CC$, that is, the space $L^2([a,b], \\CC)$. It is a vector space, and we can define an inner product on it: $\\braket{x(t)|y(t)} := \\int_a^b x^*(t)y(t)dt$. Now let's consider what the inner product of two exponentials looks like: $\\braket{e^{2\\pi if_1 t}| e^{2\\pi if_2 t}} = \\int_a^b e^{2\\pi i(f_2-f_1)t}dt = [\\frac{1}{2\\pi i(f_2-f_1)}e^{2\\pi i(f_2-f_1)t}]_a^b$. In order for this to be zero, we want $e^{2\\pi i(f_2-f_1)b} = e^{2\\pi i(f_2-f_1)a}$, that is, $2\\pi (f_2-f_1)(b-a) = 2\\pi n$ for some integer $n$. So a necessary and sufficient condition for two exponentials to be orthogonal is that $f_2-f_1 = n/(b-a)$ for some integer $n$.

This leads us to consider the set of functions $\\{ \\frac{1}{\\sqrt{b-a}} e^{2\\pi ifx} \\}_{f=n/(b-a)}$ for $n\\in\\ZZ$. It's definitely an orthonormal set of functions. Is it complete? You'll have to take my word for it, but it is in fact complete. (If the idea of trusting a complete stranger from who knows where with a result as important as this bothers you, you should take a look at a book like "Mathematics of Classical and Quantum Physics" by Byron and Fuller, which is filled with wonderful proofs of these kinds of things).

We have now decomposed any function on the interval into the sum of complex exponentials! This lets us do interesting things....

#### Physics aside: Guitar strings
Consider a string of length $L$ stretched tight between two fixed points. Suppose its mass per unit length is $\\lambda$, and it is subject to a tension of $T$.

![](/assets/tension.svg)

From the diagram above, we see that the y-component of force on the tiny bit of string is $-T\\frac{dy}{dx} + T(\\frac{dy}{dx} + \\frac{d^2y}{dx^2}dx) = T\\frac{d^2y}{dx^2}dx$. Hence, by Newton's second law $F=ma$, we have $T\\frac{d^2y}{dx^2}dx = (\\lambda dx) (\\frac{d^2y}{dt^2})$, which gives us a "wave equation": $\\frac{d^2y}{dx^2} = \\frac{1}{T/\\lambda} \\frac{d^2y}{dt^2}$.

(Aside: for comparison, here are some other wave equations:
1. Sound waves: $\\nabla^2 p = \\frac{1}{v_s^2}\\frac{\\partial^2 p}{\\partial t^2}$ where $v_s=\\sqrt{B/\\rho_0}$ is the speed of sound, $B$ is the bulk modulus of the fluid, and $\\rho_0$ is the average density.
2. EM waves (light): $\\nabla^2 \\bf{E} = \\frac{1}{c^2}\\frac{\\partial^2 \\bf{E}}{\\partial t^2}$, where $c=1/\\sqrt{\\epsilon_0\\mu_0}$ is the speed of light.
3. Matter waves (like free electrons) in quantum theory: $\\nabla^2 \\psi = \\frac{2m}{\\hbar} (-i\\frac{\\partial \\psi}{\\partial t})$. (This looks a bit different from the others; it expresses the fact that matter waves do not have a constant phase-velocity, which leads to dispersion of different frequency components. But keep in mind the $i$ is essentially a derivative, so it has sorta the same form).
4. Gravitational waves: I have no idea, but in the limit of small distortions I'm guessing it has the same form as the others.

These equations stand in contrast with the heat and diffusion equations, in which the time derivative is only first order. Hence heat and diffusion do not exhibit wave-like behavior.)

With foresight, I will define $v=\\sqrt{T/\\lambda}$.

To solve this equation, let's start by only looking for solutions of the form $y(x,t) = f(x) g(t)$. This seems a ridiculously strict assumption at first, but we will see soon that it actually yields very general solutions.

Plugging it in, we get $\\delsq{f}{x} g = \\frac{1}{v^2}f\\delsq{g}{t}$ so $\\frac{1}{f}\\delsq{f}{x} = \\frac{1}{v^2}\\frac{1}{g}\\delsq{g}{t}$. Now here is a subtle point. The right side of the equation depends only on $t$, while the right side of the equation depends only on $x$. This means that both sides are equal to some constant. If it were not so, then I could change the value of the left hand side by moving $t$ and keeping $x$ constant; because of the equality, the right side will also be changed; but this is a contradiction, because the right side depended only on $x$.

Let's call the constant $-k^2$, and let's define $\\omega := vk$. Now we have essentially decomposed our PDE into two ODEs, which are simple to solve.

$$
\\begin{align}
\\frac{1}{f}\\delsq{f}{x} &= -k^2 \\\\
\\frac{1}{g}\\delsq{g}{t} &= -\\omega^2
\\end{align}
$$

We end up with the general solutions $f(x) = Ae^{ikx} + Be^{-ikx}$ and $g(t) = Ce^{i\\omega t} + De^{-i\\omega t}$. The boundaries conditions specify that both ends of the string are fixed at zero, hence $A+B=0$ and $Ae^{ikL}+Be^{-ikL}=0$. Substituting the first equation into the second gives $A(e^{ikL}-e^{-ikL})=0 \\Rightarrow e^{ikL}=e^{-ikL}\\Rightarrow e^{2ikL}=1$, which implies $k=n\\pi/L$ for some integer $n$. We get as our final solution: $y(x,t) = \\sin(kx)e^{i\\omega t} + \\text{c.c.}$. The c.c. stands for complex-conjugate, and it is just there to make the solution real. Some people will set in place a convention that the final solution is complex, and it's implicit that the physical solution is obtained by taking the real part. This is totally fine for linear systems, and is commonly used in discussing phasors in electrical circuits, however it gets a bit confusing when non-linearity comes into play, because $\\Re(ab) \\neq \\Re(a)\\Re(b)$.

The solutions we have found are called the standing waves, because they oscillate in place. They are the first, second, third, ... harmonics, all multiples of the fundamental harmonic. They seem like special cases, until you realize that, like the Fourier series, we can decompose ANY function into a sum of these standing waves.

Normalizing the position equations gives us a family: $f(x) = \\sqrt{2} \\sin(k_n x)$ where $k_n = n\\pi /L$. These equations form a complete orthonormal basis because they are the eigenfunctions of the Hermitian operator $\\delsq{}{x}$ acting on the subspace of functions which are fixed to zero at $x=0,L$. (TODO: this fact relies on the spectral theorem. Should I include it here?)

Now for some numbers. I'm pretty clueless when it comes to music, so I'll get by with the best info I have. Looking at https://tension.stringjoy.com/, an accoustic guitar with phosphor bronze strings typically has a 25.5 inch 0.012 string (whatever 0.012 means). That string will give a pitch of E4 (329.63 Hz) when a tension of 25.7 pounds is applied. Assuming that E4 is the fundamental frequency of the note when played, we have $v=\\frac{\\omega_1}{k_1} = \\frac{2\\pi (329.63 \\mathrm{Hz})}{\\pi / (25.5 \\mathrm{in})} = 427 \\mathrm{m/s}$. Incidentally this implies that the linear mass density is $\\lambda = T/v^2 = 0.627 \\mathrm{g/m}$.

Now let's say we pluck the string at $x=0.3L$, so it looks like a triangle at $t=0$. We can decompose our initial shape $f(x)$ into a sum of our harmonics. We can attach an appropriate $e^{i\\omega t}$ to each of these harmonics to get our full solution $y(x,t)$. Here are some plots:

![](/assets/fixed-bc-power-harmonics.svg)

![](/assets/fixed-bc-full-approx.svg)

![](/assets/fixed-bc-close-approx.svg)

{{< rawhtml >}}

<video width=100% controls autoplay>
    <source src="/assets/fixed-wave.mp4" type="video/mp4">
    Your browser does not support the video tag.  
</video>

{{< /rawhtml >}}

{{< details "Periodic boundary conditions (TODO)" >}}
The problem of a guitar string fixed at its ends actually didn't involve Fourier series at all. Indeed, the Fourier series on that interval isn't suited for the desired boundary conditions, because they are not zero at $x=0,L$. Our fundamental harmonic $\\pi/L$ had a wavenumber twice that of the fundamental wavenumber of the Fourier series. However, Fourier series *are* suitable when we desire periodic boundary conditions, that is, the two ends of string are joined together so there is no more boundary. I imagine this is hard to do with strings while maintaining tension, so we can use instead a metal ring. (This sort of boundary condition is also useful when you have an endless periodic lattice of atomic nuclei, all in perfect position inside a perfect crystal, with no boundary for millimeters on end, while their free electrons zip across the interior.) In this case, the boundary conditions require $A+B = Ae^{ikL} + Be^{-ikL}$ and $ikA - ikB = ikAe^{ikL}-ikBe^{-ikL}$. Adding these two equation gives $1 = e^{ikL}$, which requires that $k = 2\\pi n/L$ for an integer $n$, which is exactly the condition of our Fourier series. But these are not enough. We need addition info on whether the wave is travelling forwards or backwards.
{{< /details >}}

The Fourier series are not the only set of complete orthonormal functions. Another beautiful example is the spherical harmonics, which form an orthonormal basis for $L^2(\\SS^2, \\CC)$ ($\\SS^2$ is the unit sphere). These show up in the orbitals of the hydrogen atom.

![](/assets/Spherical_Harmonics.png)

*Visual representations of the first few real spherical harmonics. Blue portions represent regions where the function is positive, and yellow portions represent where it is negative. The distance of the surface from the origin indicates the absolute value of the function in the angular direction. Figure and caption from Wikipedia: Spherical Harmonics*


## Fourier transforms

Fourier series can be used to model all sorts of interesting boundary problems, from the electromagnetic field in your microwave (metal walls impose zero electric field boundary conditions), to the variation of sound in an accoustic concert hall, to the behavior of an electron trapped in an infinite potential well.

They run into a limitation, however. Time is essentially infinite in both directions. Space is also infinite in all directions. How do we model light waves from distant stars, which appear to have taken millions years to reach us?

We need to consider the limit as $[a,b]$ stretches out into an infinite interval.

Let $a\\rightarrow -\\infty$ and $b\\rightarrow +\\infty$. Our inner product becomes $\\braket{x(t),y(t)} \\approx \\int_{-\\infty}^{+\\infty}x^*(t)y(t)dt$. We can still write an arbitrary function $x(t)$ as
$$
\\begin{align}
x(t) &= \\sum_{n=-\\infty}^\\infty  \\braket{\\frac{1}{\\sqrt{b-a}}e^{2\\pi i\\frac{n}{b-a}t'} | x(t')} \\frac{1}{\\sqrt{b-a}}e^{2\\pi i\\frac{n}{b-a}t} \\\\
&\\approx \\int_{f=-\\infty}^\\infty \\braket{e^{2\\pi ift'} | x(t')} e^{2\\pi ift} df \\\\
&= \\int_{f=-\\infty}^\\infty \\left(\\int_{t'=-\\infty}^\\infty e^{-2\\pi ift'} x(t') dt'\\right) e^{2\\pi ift} df
\\end{align}
$$

There we have it! We have derived our Fourier transform:
$$
\\boxed{
    \\begin{aligned}
    &X(f) = \\int_{-\\infty}^\\infty e^{-2\\pi ift} x(t)dt \\\\
    &x(t) = \\int_{-\\infty}^\\infty e^{2\\pi ift} X(f)df
    \\end{aligned}
}
$$

These equations are remarkably symmetric. In fact, it turns out that $\\mathrm{FT}[x(t)] = -\\mathrm{IFT}[x(-t)]$.

This has its uses. For example, one question that the Fourier transforms immediately raise is: what is the Fourier transform of $e^{2\\pi if_0t}$? If we try plugging it in, it yields a strange integral. For any $f\\neq f_0$, that part of the integral should be zero, since a complex exponential is, over the whole real line, zero on average. However, when $f=f_0$, the integral is infinity. This kinda makes sense since a pure tone should yield an infinitely sharp response in the frequency domain.

But we can be more precise. What is the Fourier transform of a dirac delta? $\\mathrm{FT}[\\delta (t-t_0)] = \\int_{t=-\\infty}^\\infty e^{-2\\pi ift} \\delta(t-t_0)dt = e^{-2\\pi ift_0}$. But from duality, this means $\\mathrm{FT}[e^{2\\pi if_0t}] = \\delta(f-f_0)$!

#### A different perspective: quantum mechanics and momentum space
I first learned about the Fourier transform in quantum mechanics, but it was presented in such a way that I did not even know it was a transform. Its description was so natural and intuitive that it just made sense. What follows is a quick debrief on quantum mechanics, but the math can be applied to surprisingly many other physical systems.

There are five postulates in quantum mechanics:

1 A physical system's state is a vector in a complex Hilbert space. (A Hilbert space is just an inner-product space whose norm makes it a complete metric space. All the familiar finite dimensional spaces are Hilbert spaces. $L^2(\\RR, \\CC)$ is also a Hilbert space).

2a A physical observable (e.g. energy, momentum, position, spin) is represented by a Hermitian linear operator on the space. (An operator $A$ is called Hermitian if $A=A^\\dagger$; that is, if its adjoint, aka its conjugate-transpose, is equal to itself.) When you make a measurement of the observable using any sort of apparatus, it is only possible to observe an eigenvalue of the Hermitian operator.

2b [not relevant to us] Measurement probability postulate.

2c [not relevant to us] Collapse postulate.

3 [not relevant to us] Time evolution.

For a particle that is free to travel along the entire 1-dimensional line, the Hilbert space we are interested in is $L^2(\\RR, \\CC)$.

There are two operators we are interested in. The position observable is given by the operator $f(x) \\mapsto xf(x)$, that is, multiplication by $x$. It is easy to check that this is linear and Hermitian. Its eigenvalues consist of the entire real line. An eigenvalue $c$ has eigenvector $\\delta(x-c)$. Because its eigenvectors form an orthonormal basis for the Hilbert space, we can express any wave $\\psi(x)$ in its basis. In this case it is trivial: $\\ket{\\psi} = \\int_x \\psi(x) \\ket{\\delta(x-c)} dx$.

The momentum observable, on the other hand, is the operator $\\frac{\\hbar}{i} (d/dx)$. Again, this is both linear and Hermitian. [TODO: prove momentum operator is Hermitian here?] An eigenvector $\\ket{p}$ with eigenvalue $p$ must be a solution of the differential equation $\\frac{\\hbar}{i} (d/dx)\\phi(x) = p\\phi(x)$. Its solution is $\\ket{p} = e^{ipx/\\hbar}$. Hence the eigenvectors of the momentum observable consist of the set of exponential functions. Note that $\\hbar = \\frac{h}{2\\pi}$ is the reduced plank constant, and $h\\approx 6.6\\times 10^{-34} \\mathrm{J}\\cdot\\mathrm{s}$ is the regular plank constant. If we choose to use so called "natural units" such that $h=1$, then the momentum $p$ corresponds exactly to the frequencies that we use in our official definition of the Fourier transform.

So we have two different sets of orthonormal bases for our Hilbert space. Getting from one to the other is a simple matter of changing bases. The operator $\\int_p \\ket{p}\\bra{p} dp$ is called a "resolution of the identity", because it is in fact the identity; it merely breaks our vector up as a linear combination of the basis. If we expand it out, what it is really saying is that
$$
\\ket{\\psi} = \\left(\\int_p \\ket{p}\\bra{p} dp\\right)\\ket{\\psi} = \\int_p e^{ipx/\\hbar} \\braket{p|\\psi}dp = \\int_p e^{ipx/\\hbar} \\left(\\int_{x'} e^{-ipx'/\\hbar}\\psi dx'\\right)dp
$$

But this is just our FT/IFT! At its heart, the Fourier transform is nothing more than a change in basis.
