import numpy as np 

def make_ellipse(num_points = 1000, a = 10,b=10, distrib = "uniform"):
    points = np.zeros((num_points,3))
    normals = np.zeros((num_points,3))
    curvature = np.zeros(num_points)
    if distrib == 'uniform':
        phi = np.linspace(0,2*np.pi,num_points+1)
    elif distrib == 'gauss':
        phi = np.random.normal(size = num_points+1)
        phi = (phi - np.min(phi))*2*np.pi/(np.max(phi) - np.min(phi))
        phi = np.sort(phi)
    elif distrib == 'unipert':
        phi = np.linspace(0,2*np.pi,num_points+1)
        dphi = (phi[1]-phi[0])/2
        phi = phi + np.random.uniform(low = 0,high=dphi,size = num_points+1)
    points[...,0] = a * np.sin(phi[:-1])
    points[...,1] = b * np.cos(phi[:-1])
    normals[...,0] = b*np.sin(phi[:-1])
    normals[...,1] = a*np.cos(phi[:-1])
    curvature = a*b/np.sqrt((a**2*(np.cos(phi[:-1]))**2 + b**2*(np.sin(phi[:-1]))**2)**3)

    
    return points, normals, curvature

def make_spiral(num_points = 1000, k = 1 , max_phi = np.pi, distrib = "uniform" ):
    points = np.zeros((num_points,3))
    normals = np.zeros((num_points,3))
    curvature = np.zeros(num_points)
    if distrib == 'uniform':
        phi = np.linspace(0.00001,max_phi,num_points)
    elif distrib == 'gauss':
        phi = np.random.normal(size = num_points)
        phi = (phi - np.min(phi))*max_phi/(np.max(phi) - np.min(phi))
        phi = np.sort(phi)
        if(phi[0] == 0):
            phi[0] = 0.00001
    elif distrib == 'unipert':
        phi = np.linspace(0.00001,max_phi,num_points)
        dphi = (phi[1]-phi[0])/2
        new_phi = phi + np.random.uniform(low = 0,high=dphi,size = num_points)
        phi[1:-1] = new_phi[1:-1]
    points[...,0] = k*phi * np.sin(phi)
    points[...,1] = k*phi * np.cos(phi)
    normals[...,0] = -k*(np.cos(phi) -phi*np.sin(phi))
    curvature = (phi**2 + 2)/(np.sqrt((phi**2 + 1)**3)*k)
    normals[...,1] = k*(np.sin(phi)+phi*np.cos(phi))
    
    return points, normals, curvature

def make_delta(num_points = 1000, a = 4, b = 12, max_phi = 2*np.pi, distrib = "uniform"):
    assert a != b
    points = np.zeros((num_points,3))
    normals = np.zeros((num_points,3))
    curvature = np.zeros(num_points)
    if distrib == 'uniform':
        phi = np.linspace(0,2*np.pi,num_points+1)
    elif distrib == 'gauss':
        phi = np.random.normal(size = num_points+1)
        phi = (phi - np.min(phi))*max_phi/(np.max(phi) - np.min(phi))
        phi = np.sort(phi)
    elif distrib == 'unipert':
        phi = np.linspace(0,2*np.pi,num_points+1)
        dphi = (phi[1]-phi[0])/2
        phi = phi + np.random.uniform(low = 0,high=dphi,size = num_points+1)
    indxs = np.where(np.isclose(np.cos(b*phi[:-1]/a), 1))
    phi[indxs]+=0.00001
    points[...,0] = (b-a)*np.cos(phi[:-1]) + a*np.cos((b-a)*phi[:-1]/a)
    points[...,1] = (b-a)*np.sin(phi[:-1]) - a*np.sin((b-a)*phi[:-1]/a)
    normals[...,0] = (b-a)*(np.cos(phi[:-1])-np.cos((b-a)*phi[:-1]/a))
    curvature = ((b-a)**2-((b-a)**3)/a) * (1 - np.cos(b*phi[:-1]/a))/(np.sqrt(2*(b-a)**2)*(1-np.cos(b*phi[:-1]/a)))**3
    normals[...,1] = (b-a)*(np.sin(phi[:-1])+np.sin((b-a)*phi[:-1]/a))
    
    return points, normals, curvature

def make_loop(num_points = 1000, t0 = -3, t1 = 6, distrib = "uniform" ):
    points = np.zeros((num_points,3))
    normals = np.zeros((num_points,3))
    curvature = np.zeros(num_points)
    if distrib == 'uniform':
        t = np.linspace(t0,t1,num_points)
    elif distrib == 'gauss':
        t = np.random.normal(size = num_points)
        t = t0 + (t - np.min(t))*(t1 - t0)/(np.max(t) - np.min(t))
        t = np.sort(t)
    elif distrib == 'unipert':
        t = np.linspace(t0,t1,num_points)
        dt = (t[1]-t[0])/2
        new_t = t + np.random.uniform(low = 0,high=dt,size = num_points)
        t[1:-1] = new_t[1:-1]
    t[np.where(t == 0)] = 0.0001
    t[np.where(t == 2)] = 2+ 0.0001
    points[...,0] = t*(2-t)
    points[...,1] = t**2 *(2-t)
    normals[...,0] = 4*t-3*t**2
    curvature = 2*(3*t**3 - 6*t+4)/(np.sqrt(9*t**4 -24*t**3 + 20*t**2 -8*t +4))**3
    normals[...,1] = - 2*(1-t)
    
    return points, normals, curvature

def make_parabola(num_points = 1000, t0 = -5 , t1 = 5, distrib = "uniform" ):
    points = np.zeros((num_points,3))
    normals = np.zeros((num_points,3))
    curvature = np.zeros(num_points)
    if distrib == 'uniform':
        t = np.linspace(t0,t1,num_points)
    elif distrib == 'gauss':
        t = np.random.normal(size = num_points)
        t = t0 + (t - np.min(t))*(t1- t0)/(np.max(t) - np.min(t))
        t = np.sort(t)
    elif distrib == 'unipert':
        t = np.linspace(t0,t1,num_points)
        dt = (t[1]-t[0])/2
        t = t + np.random.uniform(low = 0,high=dt,size = num_points)
    t[np.where(t == 0)] = 0.0001
    points[...,0] = t
    points[...,1] = t**2 
    normals[...,0] = 2*t
    curvature = 2/np.sqrt(4*t**2 +1)**3
    normals[...,1] = -1
    
    return points, normals, curvature
