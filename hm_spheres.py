class Ball():
  def __init__(self, mass, radius, name, pos=[], vel=[]):
    self.mass = mass
    self.radius = radius
    self.name = name
    self.pos = np.array(pos, np.float64)
    self.vel = np.array(vel, np.float64)
    self.time = 0.00000000

  def __repr__(self):
    return str(self.name)+" "+str(self.pos)[1:-1]+" "+str(self.vel)[1:-1]

  def updatepos(self, t): 
    # dt = (t - self.time)
    dt = float(Decimal(t) - Decimal(self.time))
    self.pos += dt * self.vel
    self.time = t
    
    
class Universe():
    def __init__(self, radius, max_col, ball_array):
    #constructor
        self.radius=radius
        self.max_col=max_col
        self.ball_array=copy.deepcopy(ball_array)
        self.time = 0.00000000
        
    def __repr__(self):
		str_out = str(self.time) + "\n"
		for i in range(len(self.ball_list)):
			str_out += str(self.ball_list[i].name)+" "+str(self.ball_list[i].pos)[1:-1] 
									+" "+str(self.ball_list[i].vel)[1:-1]
			if i != len(self.ball_list) - 1:
				str_out += "\n"
		return str_out   
    
    '''    
    def isCollision(self, ball_a, ball_b):
		p1, p2 = ball_a.pos, ball_b.pos
		v1, v2 = ball_a.vel, ball_b.vel
		ans = np.dot(p1-p2,v1-v2)
		if(ans<0):
			return True
		else:
			return False
     '''
     
    def Ucollision(self, ball1):    
        magp=np.linalg.norm(ball1.pos)
        upos=np.array(ball1.pos)*(1/magp)
        
        norm_v = np.dot(ball1.vel,upos)*np.array(upos)
        tan_v = np.subtract(ball1.vel,norm_v)
        ball1.vel = np.subtract(tan_v,norm_v)
        ball1.bounce+=1
        
        
    def collision_time():
    #collision time
    
    def collisions():
    #collisions    
        
    
    def update_pos(self, t):
    #update all ball positions
    for ball in self.ball_array:
		ball.updatepos(t)
    
        
    def update_vel(self, ball1, ball2):
    #update ball velocity
    p1, p2= ball1.pos, ball2.pos
    v1, v2= ball1.vel, ball2.vel
    m1, m2 = ball1.mass, ball2.mass
    
    vp=  np.subtract(v1,v2)
    v_p= np.subtract(v2-v1)
    
    rp=  np.subtract(p1-p2)
    r_p= np.subtract(p2-p1)
    
    mp=  np.subtract(m1-m2)
    m_p= p.subtract(m2-m1)
    
    ball1.vel= v1 - (2.*m2/(m1+m2)) * (np.dot(vp,rp)/np.dot(rp,rp))*rp 
    ball2.vel= v2 - (2.*m1/(m1+m2)) * (np.dot(v_p,r_p)/np.dot(r_p,r_p))*r_p 
    
    ball1.bounce+=1
    ball2.bounce+=1
    
    
    
    def update_time():
    #update time
