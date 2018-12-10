from Reconstruct_States import minmod


x = -1
y = 2
z = 3

# should be 1
print("Minmod: " , minmod(x,y,z))

ci = 0.5
ci1 = 1
ci2 = 3
theta = 1.
# should be 0.75
print(1-0.5*minmod(theta*(ci1-ci),0.5*(ci2-ci),theta*(ci2-ci1)))
