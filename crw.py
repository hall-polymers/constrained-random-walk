# !/usr/bin/env pypy

### edited J.Brown 20140404 modified Kevin Shen 20170711/20200803

import sys
from math import pi, sin, cos, sqrt, floor
from random import random, randint, seed
from scipy import optimize
import numpy as np
import argparse

def crw(M, N, L, constraint, outputtraj, intrvl, fA, x0=-1, y0=-1, z0=-1):
	# x0, y0, z0 - the inital location of the RW, a negative number will choose random location between 0 and L FIXME: does not check that these choices are reasonable, so an infinite loop can happen if you force the system to start inside a particle
	# set the starting locations for the RWs
	if x0 < 0:
		xstart = lambda: L*random()
	else:
		xstart = lambda: x0
	if y0 < 0:
		ystart = lambda: L*random()
	else:
		ystart = lambda: y0
	if z0 < 0:
		zstart = lambda: L*random()
	else:
		zstart = lambda: z0

	dumpArr = [0]* int(N//intrvl)

	# for gyroid surface function
	rescale = 2*pi/L
	gy_crit = (1-fA)/0.067
	gy1 = lambda x,y,z : 10.0*(cos(x*rescale) * sin(y*rescale) + cos(y*rescale) * sin(z*rescale) + cos(z*rescale) * sin(x*rescale)) - 0.5*(cos(2*x*rescale)*cos(2*y*rescale) + cos(2*y*rescale)*cos(2*z*rescale) + cos(2*z*rescale)*cos(2*x*rescale)) - gy_crit
	gy2 = lambda x,y,z : -10.0*(cos(x*rescale) * sin(y*rescale) + cos(y*rescale) * sin(z*rescale) + cos(z*rescale) * sin(x*rescale)) - 0.5*(cos(2*x*rescale)*cos(2*y*rescale) + cos(2*y*rescale)*cos(2*z*rescale) + cos(2*z*rescale)*cos(2*x*rescale)) - gy_crit

	def gyroid(d, in_gy1, r0, r1, r2, dx, dy, dz):
		x = r0+dx*d
		y = r1+dy*d
		z = r2+dz*d
		if in_gy1:
			return round(gy1(x, y ,z), 10)
		return round(gy2(x, y ,z), 10)

	r = [0.0, 0.0, 0.0]
	D = 0

	if outputtraj:
		filename = "crw.lammpstrj"
		output = open(filename, 'w')
		output.write("ITEM: TIMESTEP\n")
		output.write("0\n")
		output.write("ITEM: NUMBER OF ATOMS\n")
		output.write("%10i\n" % (M*(N+1)))
		output.write("ITEM: BOX BOUNDS pp pp pp\n")
		output.write("%10i %10i\n" % (0, L))
		output.write("%10i %10i\n" % (0, L))
		output.write("%10i %10i\n" % (0, L))
		output.write("ITEM: ATOMS id mol type x y z ix iy iz\n")
	
	if constraint in ["none", "lamellae"]: # bulk (no constraint) or reflective walls of lamellae constraint
		# loop over all the chains
		for m in xrange(M):
			# place the starting location 
			ir = [0, 0, 0]
			r[0] = xstart()
			r[1] = ystart()
			r[2] = zstart()

			if outputtraj:
				output.write("%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n" % (m*(N+1)+1, m+1, 1, r[0], r[1], r[2], 0, 0, 0))
		
			r0 = r[:]

			# loop over all the links in the chain
			for n in xrange(1, N+1):
				# choose a random point on a unit sphere
				dz = 2*random()-1.0
				theta = 2*pi*random()
				ds = sqrt(1-dz*dz)
				dx = ds*cos(theta)
				dy = ds*sin(theta)
		
				r[0] += dx
				r[1] += dy
				r[2] += dz

				if constraint == "lamellae": # walls at z=0 and z=L
					#if we ended up past a wall, enforce the reflective BC
					if r[2] < 0 or r[2] > L:
						done = False
						while not done:
							if r[2] < 0:
								r[2] = -r[2]
							elif r[2] > L:
								r[2] = 2*L - r[2] 
							else:
								done = True

				# enforce periodic BCs
				ir[0] += floor(r[0]/L)
				r[0] = L*((r[0]/L)%1)
				ir[1] += floor(r[1]/L)
				r[1] = L*((r[1]/L)%1)
				ir[2] += floor(r[2]/L)
				r[2] = L*((r[2]/L)%1)
				#print ir, r

				if outputtraj:
					output.write("%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n" % (m*(N+1)+1 + n, m+1, 1, r[0], r[1], r[2], int(ir[0]), int(ir[1]), int(ir[2])))

				if n % intrvl == 0:
					dumpArr[n//intrvl-1] += ((L*ir[0]+r[0]-r0[0])**2 + (L*ir[1]+r[1]-r0[1])**2 + (L*ir[2]+r[2]-r0[2])**2)/M
	
	elif constraint in ["cylinder"]: # cylindrical reflective wall
		# loop over all the chains
		for m in xrange(M):
			# place the starting location, as long as it's not outside the cylinder
			ir = [0, 0, 0]
			r[0] = xstart()
			r[1] = ystart()
			r[2] = zstart()
			while (r[0] - L/2)**2 + (r[1] - L/2)**2 > (L/2)**2:
				r[0] = xstart()
				r[1] = ystart()
				r[2] = zstart()

			if outputtraj:
				output.write("%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n" % (m*(N+1)+1, m+1, 1, r[0], r[1], r[2], 0, 0, 0))
		
			r0 = r[:]
			
			# loop over all the links in the chain
			for n in xrange(1, N+1):
				# choose a random point on a unit sphere
				dz = 2*random()-1.0
				theta = 2*pi*random()
				ds = sqrt(1-dz*dz)
				dx = ds*cos(theta)
				dy = ds*sin(theta)
				dx_0 = ds*cos(theta)
				dy_0 = ds*sin(theta)
		
				# check for reflections as long as there is distance remaining to be traveled
				dist_remaining = ds
				reflection = True
				# loop over reflections
				# i = 0
				while reflection:
					reflection = False

					xc = L/2
					yc = L/2
					# compute the distance to the cylinder center, allowing quick disqualification
					centerdistsq = (r[0]-xc)**2 + (r[1]-yc)**2

					# make xy vector to be unit vector for calculation
					rsqrt_xy = sqrt(dx*dx+dy*dy)
					dx = dx/rsqrt_xy
					dy = dy/rsqrt_xy
					
					# if close enough from the center initially, we can ignore it, otherwise have to check reflections
					if centerdistsq < ((L/2)-1)**2:
						# default to infinite distance
						d = float("inf")
						continue

					# we're close enough to the cylindrical surface
					
					# find the intersection with the sphere, if there is one, using the quadratic formula
					# see: http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
					a = dx**2 + dy**2
					b = 2*(dx*(r[0]-xc) + dy*(r[1]-yc))
					c = (r[0]-xc)**2 + (r[1]-yc)**2 - (L/2)**2
					radicand = b**2 - 4*a*c
	
					# if the term under the radical is negative, no intersection, otherwise test further for a reflection
					if radicand < 0:
						# default to infinite distance
						print >> sys.stderr, "Error: A bead has been out of the cylinder, should be impossible."
						return

					# need to find the point of reflection and move from there
					d1 = (-b - sqrt(radicand))/(2*a)
					d2 = (-b + sqrt(radicand))/(2*a)
	
					# check if the intersection is outside of our range, if it is we can just move on
					# this is for the case where the floating number makes it slightly outside of the cylinder
					if d1 > 0 and d2 > 0:
						if d1 > d2:
							d = d1
						else:
							d = d2
					# normal cases
					elif d1 > 0 and d1 < dist_remaining:
						d = d1
					elif d2 > 0 and d2 < dist_remaining:
						d = d2
					else:
						# default to infinite distance
						d = float("inf")

					if d < dist_remaining:
						# found a reflection!
						dist_remaining = dist_remaining - d
						
						r[0] += d*dx
						r[1] += d*dy

						# still have to move l-d in the reflected direction (simplified by the fact that the current r is on the sphere)
						# see http://en.wikipedia.org/wiki/Reflection_%28mathematics%29
						drdotr = dx*(r[0]-xc) + dy*(r[1]-yc)
						dx = dx - 2*drdotr*(r[0]-xc)/((L/2)**2)
						dy = dy - 2*drdotr*(r[1]-yc)/((L/2)**2)
		
						# the second movement we don't explicitly perform, instead we just change dx, dy, and dz appropriately, and restart the loop
						reflection = True

	
				if dist_remaining == ds:
					r[0] += dx_0
					r[1] += dy_0
					r[2] += dz
				else:
					r[0] += dist_remaining*dx
					r[1] += dist_remaining*dy
					r[2] += dz

				# the vector that we check against for the nonreversing restriction is the last direction we moved (i.e. the reflected direction if there was a reflection)
				dx_prev = dx
				dy_prev = dy
				dz_prev = dz
	
				# enforce periodic BCs
				ir[0] += floor(r[0]/L)
				r[0] = L*((r[0]/L)%1)
				ir[1] += floor(r[1]/L)
				r[1] = L*((r[1]/L)%1)
				ir[2] += floor(r[2]/L)
				r[2] = L*((r[2]/L)%1)
				
				if outputtraj:
					output.write("%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n" % (m*(N+1)+1 + n, m+1, 1, r[0], r[1], r[2], int(ir[0]), int(ir[1]), int(ir[2])))

				if ir[0] != 0 or ir[1] != 0:
					# image flag of x and y should not be other than 0
					print >> sys.stderr, "Image Flag Warning: A bead has been out of the cylinder, should be impossible."
					return 
	
				if n % intrvl == 0:
					dumpArr[n//intrvl-1] += ((L*ir[0]+r[0]-r0[0])**2 + (L*ir[1]+r[1]-r0[1])**2 + (L*ir[2]+r[2]-r0[2])**2)/M
	
	elif constraint in ["gyroid"]: # gyroid reflective wall
		# loop over all the chains
		for m in xrange(M):
			# place the starting location, as long as it's not outside the cylinder
			ir = [0, 0, 0]
			r[0] = xstart()
			r[1] = ystart()
			r[2] = zstart()
			while max(gy1(r[0], r[1], r[2]), gy2(r[0], r[1], r[2])) < 0:
				r[0] = xstart()
				r[1] = ystart()
				r[2] = zstart()

			if outputtraj:
				output.write("%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n" % (m*(N+1)+1, m+1, 1, r[0], r[1], r[2], 0, 0, 0))

			r0 = r[:]

			# determine which network it is in
			in_gy1 = False
			if gy1(r[0], r[1], r[2]) > 0:
				in_gy1 = True
			
			# loop over all the links in the chain
			for n in xrange(1, N+1):
				# choose a random point on a unit sphere
				dz = 2*random()-1.0
				theta = 2*pi*random()
				ds = sqrt(1-dz*dz)
				dx = ds*cos(theta)
				dy = ds*sin(theta)

				# check for reflections as long as there is distance remaining to be traveled
				dist_remaining = 1.0
				reflection = True
				while reflection:
					reflection = False

					if gyroid(0, in_gy1, r[0], r[1], r[2], dx, dy, dz) < 0:
						print >> sys.stderr, "Error: A bead has been out of the gyroid, should be impossible."
						return 

					elif gyroid(0, in_gy1, r[0], r[1], r[2], dx, dy, dz) > 0.1 and gyroid(dist_remaining, in_gy1, r[0], r[1], r[2], dx, dy, dz) > 0.1:
						d = float("inf")
						continue

					elif min(optimize.fmin_l_bfgs_b(gyroid, dist_remaining, args=(in_gy1, r[0], r[1], r[2], dx, dy, dz), bounds=[(0,dist_remaining)], approx_grad=True)[1], optimize.fmin_l_bfgs_b(gyroid, 0, args=(in_gy1, r[0], r[1], r[2], dx, dy, dz), bounds=[(0,dist_remaining)], approx_grad=True)[1]) < 0:
						init_position = 0
						final_position = dist_remaining
						
						if gyroid(0, in_gy1, r[0], r[1], r[2], dx, dy, dz) <= 0 and gyroid(dist_remaining, in_gy1, r[0], r[1], r[2], dx, dy, dz) < 0:
							init_position = optimize.minimize_scalar(lambda d: -gyroid(d, in_gy1, r[0], r[1], r[2], dx, dy, dz), bounds=(0, dist_remaining), method='bounded').x

						elif gyroid(dist_remaining, in_gy1, r[0], r[1], r[2], dx, dy, dz) >= 0:
							final_position = optimize.minimize_scalar(gyroid, bounds=(0, dist_remaining), method='bounded', args=(in_gy1, r[0], r[1], r[2], dx, dy, dz)).x

						d = optimize.brentq(gyroid, init_position, final_position, args=(in_gy1, r[0], r[1], r[2], dx, dy, dz))

						if d < dist_remaining:
							# found a reflection!
							dist_remaining = dist_remaining - d
							r[0] += d*dx
							r[1] += d*dy
							r[2] += d*dz

							# still have to move l-d in the reflected direction 
							gradient_x = (-10*sin(r[0]*rescale)*sin(r[1]*rescale)+10*cos(r[2]*rescale)*cos(r[0]*rescale)+sin(2*r[0]*rescale)*cos(2*r[1]*rescale)+cos(2*r[2]*rescale)*sin(2*r[0]*rescale))*rescale 
							gradient_y = (-10*sin(r[1]*rescale)*sin(r[2]*rescale)+10*cos(r[0]*rescale)*cos(r[1]*rescale)+sin(2*r[1]*rescale)*cos(2*r[2]*rescale)+cos(2*r[0]*rescale)*sin(2*r[1]*rescale))*rescale 
	 						gradient_z = (-10*sin(r[2]*rescale)*sin(r[0]*rescale)+10*cos(r[1]*rescale)*cos(r[2]*rescale)+sin(2*r[2]*rescale)*cos(2*r[0]*rescale)+cos(2*r[1]*rescale)*sin(2*r[2]*rescale))*rescale 
	 						if not in_gy1:
	 							gradient_x = (10*sin(r[0]*rescale)*sin(r[1]*rescale)-10*cos(r[2]*rescale)*cos(r[0]*rescale)+sin(2*r[0]*rescale)*cos(2*r[1]*rescale)+cos(2*r[2]*rescale)*sin(2*r[0]*rescale))*rescale 
								gradient_y = (10*sin(r[1]*rescale)*sin(r[2]*rescale)-10*cos(r[0]*rescale)*cos(r[1]*rescale)+sin(2*r[1]*rescale)*cos(2*r[2]*rescale)+cos(2*r[0]*rescale)*sin(2*r[1]*rescale))*rescale 
		 						gradient_z = (10*sin(r[2]*rescale)*sin(r[0]*rescale)-10*cos(r[1]*rescale)*cos(r[2]*rescale)+sin(2*r[2]*rescale)*cos(2*r[0]*rescale)+cos(2*r[1]*rescale)*sin(2*r[2]*rescale))*rescale
							vdota_over_adota = (dx*gradient_x+dy*gradient_y+dz*gradient_z)/(gradient_x*gradient_x+gradient_y*gradient_y+gradient_z*gradient_z)
							dx = dx - 2*vdota_over_adota*gradient_x
							dy = dy - 2*vdota_over_adota*gradient_y
							dz = dz - 2*vdota_over_adota*gradient_z

							# the second movement we don't explicitly perform, instead we just change dx, dy, and dz appropriately, and restart the loop
							reflection = True

				# move int he dx, dy, dz direction
				r[0] += dist_remaining*dx
				r[1] += dist_remaining*dy
				r[2] += dist_remaining*dz

				# the vector that we check against for the nonreversing restriction is the last direction we moved (i.e. the reflected direction if there was a reflection)
				dx_prev = dx
				dy_prev = dy
				dz_prev = dz

				# enforce periodic BCs
				ir[0] += floor(r[0]/L)
				r[0] = L*((r[0]/L)%1)
				ir[1] += floor(r[1]/L)
				r[1] = L*((r[1]/L)%1)
				ir[2] += floor(r[2]/L)
				r[2] = L*((r[2]/L)%1)

				if outputtraj:
					output.write("%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n" % (m*(N+1)+1 + n, m+1, 1, r[0], r[1], r[2], int(ir[0]), int(ir[1]), int(ir[2])))
				
				if n % intrvl == 0:
					print n, n//intrvl-1
					dumpArr[n//intrvl-1] += ((L*ir[0]+r[0]-r0[0])**2 + (L*ir[1]+r[1]-r0[1])**2 + (L*ir[2]+r[2]-r0[2])**2)/M
	
	return dumpArr

def main(args):
	# Simulate random walks
	repeat = args.ntrial
	totResult = []
	for i in range(repeat):
		result = crw(args.nparticle, args.nstep, args.lbox, args.constraint, args.outputtraj, args.dumpintrvl, args.fA)
		totResult.append(result)

	# Tally final mean square displacement results
	avg = np.mean(totResult, axis=0)
	std = np.std(totResult, axis=0)

	# Print title and results
	print "   step      msd-avg      msd-std"
	print "%7i %12.4f %12.4f" % (0, 0, 0)
	for i in range(args.nstep//args.dumpintrvl):
		print "%7i %12.4f %12.4f" % ((i+1)*args.dumpintrvl, avg[i], std[i])

def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--ntrial', type=int, default=10, help='Number of trials')
	parser.add_argument('-p', '--nparticle', type=int, default=1, help='Number of particles in each trial')
	parser.add_argument('-s', '--nstep', type=int, default=50000, help='Number of steps for each paricle in each trial')
	parser.add_argument('-l', '--lbox', type=int, default=20, help='Box dimension (cubic box)')
	parser.add_argument('-c', '--constraint', type=str, default='none', help='Constraint type')
	parser.add_argument('-o', '--outputtraj', action='store_true', default=False, help='Print trajectory file')
	parser.add_argument('-i', '--dumpintrvl', type=int, default=2500, help='MSD result interval')
	parser.add_argument('-f', '--fA', type=float, default=0.35, help='Volume fraction of conducting (A) domain, only important to gyroid phase')
	args = parser.parse_args()
	return args

if __name__ == "__main__":
	args = parse_arguments()
	sys.exit(main(args))
