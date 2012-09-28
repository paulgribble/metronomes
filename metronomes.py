# ipython --pylab

# Paul Gribble 2012

from scipy.integrate import odeint

def metronomes(state,t,params):
	m1 = params['m1']
	m2 = params['m2']
	m3 = params['m3']
	m4 = params['m4']
	m5 = params['m5']
	l1 = params['l1']
	l2 = params['l2']
	l3 = params['l3']
	l4 = params['l4']
	l5 = params['l5']
	M  = params['M']
	g  = params['g']
	Mk = params['Mk']
	Md = params['Md']
	adamp = params['adamp']
	aescp = params['aescp']
	a1 = state[0]   # angle of pendulum 1
	a2 = state[1]   # angle of pendulum 2
	a3 = state[2]   # angle of pendulum 3
	a4 = state[3]   # angle of pendulum 4
	a5 = state[4]   # angle of pendulum 5
	a1d = state[5]  # vel of pendulum 1
	a2d = state[6]  # vel of pendulum 2
	a3d = state[7]  # vel of pendulum 3
	a4d = state[8]  # vel of pendulum 4
	a5d = state[9]  # vel of pendulum 5
	x = state[10]   # position of cart
	xd = state[11]  # vel of cart
	tmp1 = (m1/M)*g*sin(a1) + (m1/M)*l1*a1d*a1d*sin(a1)
	tmp2 = (m2/M)*g*sin(a2) + (m2/M)*l2*a2d*a2d*sin(a2)
	tmp3 = (m3/M)*g*sin(a3) + (m3/M)*l3*a3d*a3d*sin(a3)
	tmp4 = (m4/M)*g*sin(a4) + (m4/M)*l4*a4d*a4d*sin(a4)
	tmp5 = (m5/M)*g*sin(a5) + (m5/M)*l5*a5d*a5d*sin(a5)
	xdd_top = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 - (Md/M)*xd - (Mk/M)*x
	xdd_bot = 1.0 - (m1/M)*cos(a1) - (m2/M)*cos(a2) - (m3/M)*cos(a3) - (m4/M)*cos(a4) - (m5/M)*cos(a5)
	xdd = xdd_top / xdd_bot
	a1dd = -(g/l1)*sin(a1) - (xdd/l1)*cos(a1) - adamp*a1d*(((a1/aescp)*(a1/aescp))-1)
	a2dd = -(g/l2)*sin(a2) - (xdd/l2)*cos(a2) - adamp*a2d*(((a2/aescp)*(a2/aescp))-1)
	a3dd = -(g/l3)*sin(a3) - (xdd/l3)*cos(a3) - adamp*a3d*(((a3/aescp)*(a3/aescp))-1)
	a4dd = -(g/l4)*sin(a4) - (xdd/l4)*cos(a4) - adamp*a4d*(((a4/aescp)*(a4/aescp))-1)
	a5dd = -(g/l5)*sin(a5) - (xdd/l5)*cos(a5) - adamp*a5d*(((a5/aescp)*(a5/aescp))-1)
	return [a1d, a2d, a3d, a4d, a5d, a1dd, a2dd, a3dd, a4dd, a5dd, xd, xdd]

# parameters
params = {'m1' : 0.02,        # mass of the bob (kg)
          'm2' : 0.02,        # mass of the bob (kg)
          'm3' : 0.02,        # mass of the bob (kg)
          'm4' : 0.02,        # mass of the bob (kg)
          'm5' : 0.02,        # mass of the bob (kg)
          'l1' : 0.25,        # length of the pendulum (m)
          'l2' : 0.25,        # length of the pendulum (m)
          'l3' : 0.25,        # length of the pendulum (m)
          'l4' : 0.25,        # length of the pendulum (m)
          'l5' : 0.25,        # length of the pendulum (m)
          'M'  : 1.00,        # mass of the cart
          'g'  : 9.81,        # graviational constant
          'Mk' : 5.00,        # stiffness coefficient for cart
          'Md' : 0.001,       # damping coefficient for cart
          'adamp' : 0.50,     # escapement damping coefficient
		  'aescp' : 30*pi/180 # escapement angle coefficient
}

# initial conditions
ar = random.uniform(-30, 30, 5)*pi/180
ard = random.uniform(-360, 360, 5)*pi/180
state0 = concatenate((ar, ard, array([0.0, 0.0])))
t = arange(0, 30, 0.01)
state = odeint(metronomes, state0, t, args=(params,))

figure(figsize=(8,10))
subplot(3,1,1)
plot(t,state[:,[0,1,2,3,4]]*180/pi)
ylabel('METRONOME ANGLES (rad)')
xlabel('TIME (sec)')
subplot(3,2,3)
plot(t[0:200],state[0:200,[0,1,2,3,4]]*180/pi)
subplot(3,2,4)
plot(t[-200:-1],state[-200:-1,[0,1,2,3,4]]*180/pi)
subplot(3,1,3)
plot(t,state[:,10])
ylabel('CART POSITION (m)')
xlabel('TIME (sec)')

# animation
def animate(state, t):
	figure()
	axis('equal')
	lw = 3
	b1x = state[0,10] + 0.0
	b1y = 0.0
	p1x = b1x + params['l1']*cos(state[0,0]-(pi/2))
	p1y = b1y + params['l1']*sin(state[0,0]-(pi/2))
	l1, = plot([b1x, p1x], [b1y, p1y], 'b-', linewidth=lw)

	b2x = state[0,10] + 0.05
	b2y = 0.0
	p2x = b2x + params['l2']*cos(state[0,1]-(pi/2))
	p2y = b2y + params['l2']*sin(state[0,1]-(pi/2))
	l2, = plot([b2x, p2x], [b2y, p2y], 'b-', linewidth=lw)

	b3x = state[0,10] + 0.10
	b3y = 0.0
	p3x = b3x + params['l3']*cos(state[0,2]-(pi/2))
	p3y = b3y + params['l3']*sin(state[0,2]-(pi/2))
	l3, = plot([b3x, p3x], [b3y, p3y], 'b-', linewidth=lw)

	b4x = state[0,10] + 0.15
	b4y = 0.0
	p4x = b4x + params['l4']*cos(state[0,3]-(pi/2))
	p4y = b4y + params['l4']*sin(state[0,3]-(pi/2))
	l4, = plot([b4x, p4x], [b4y, p4y], 'b-', linewidth=lw)

	b5x = state[0,10] + 0.20
	b5y = 0.0
	p5x = b5x + params['l5']*cos(state[0,4]-(pi/2))
	p5y = b5y + params['l5']*sin(state[0,4]-(pi/2))
	l5, = plot([b5x, p5x], [b5y, p5y], 'b-', linewidth=lw)

	c, = plot([b1x, b5x],[0,0],'k-',linewidth=lw)
	xlim(min(state[:,10])-0.3, max(state[:,10])+.5)
	tt = title("%6.2f sec" % t[0])
	for i in range(0,size(t),4):
		b1x = state[i,10] + 0.00
		b2x = state[i,10] + 0.05
		b3x = state[i,10] + 0.10
		b4x = state[i,10] + 0.15
		b5x = state[i,10] + 0.20
		p1x = b1x + params['l1']*cos(state[i,0]-(pi/2))
		p2x = b2x + params['l2']*cos(state[i,1]-(pi/2))
		p3x = b3x + params['l3']*cos(state[i,2]-(pi/2))
		p4x = b4x + params['l4']*cos(state[i,3]-(pi/2))
		p5x = b5x + params['l5']*cos(state[i,4]-(pi/2))
		p1y = b1y + params['l1']*sin(state[i,0]-(pi/2))
		p2y = b2y + params['l2']*sin(state[i,1]-(pi/2))
		p3y = b3y + params['l3']*sin(state[i,2]-(pi/2))
		p4y = b4y + params['l4']*sin(state[i,3]-(pi/2))
		p5y = b5y + params['l5']*sin(state[i,4]-(pi/2))
		l1.set_xdata([b1x, p1x])
		l2.set_xdata([b2x, p2x])
		l3.set_xdata([b3x, p3x])
		l4.set_xdata([b4x, p4x])
		l5.set_xdata([b5x, p5x])
		l1.set_ydata([b1y, p1y])
		l2.set_ydata([b2y, p2y])
		l3.set_ydata([b3y, p3y])
		l4.set_ydata([b4y, p4y])
		l5.set_ydata([b5y, p5y])
		c.set_xdata([b1x, b5x])
		tt.set_text("%6.2f sec" % t[i])
		draw()

# end

