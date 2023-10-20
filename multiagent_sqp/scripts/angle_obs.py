import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class Obstacle:
    x0: float
    y0: float
    x1: float
    y1: float
    def swap(self):
        self.x0, self.y0, self.x1, self.y1 = self.x1, self.y1, self.x0, self.y0

@dataclass
class PathSegment:
    x0: float
    y0: float
    x1: float
    y1: float
    def swap(self):
        self.x0, self.y0, self.x1, self.y1 = self.x1, self.y1, self.x0, self.y0


def isLeft(obs, path):
    return (path.x1 - path.x0) * (obs.y0 - path.y0) - (path.y1 - path.y0) * (obs.x0 - path.x0) > 0

def get_angle(obs, path, is_left):
    dx_obs = obs.x1 - obs.x0
    dy_obs = obs.y1 - obs.y0

    dx_path = path.x1 - path.x0
    dy_path = path.y1 - path.y0
    print(obs.x0, obs.y0, obs.x1, obs.y1)
    print(path.x0, path.y0, path.x1, path.y1)
    cos_term = (dx_obs * dx_path + dy_obs * dy_path)/np.sqrt((dx_obs**2 + dy_obs**2) * (dx_path**2 + dy_path**2))
    ang = abs(np.arccos(cos_term))
    cross = dx_path * dy_obs - dy_path * dx_obs
    direction = 0
    negative = False
    if cross < 0:
        ang = -ang
        negative = True

    if is_left and ang < 0:
        direction = 1  # this means the x1 point is closer to the path than x0
    elif is_left and ang > 0:
        direction = -1
    elif not is_left and ang < 0:
        direction = -1
    elif not is_left and ang > 0:
        direction = 1

    if ang < 0 and is_left:
        ang = ang + np.pi
    elif ang > 0 and not is_left:
        ang = ang - np.pi

    # if ang < 0 and is_left:
    #     ang = ang + np.pi
    # elif ang > 0 and not is_left:
    #     ang = ang - np.pi
    ang = abs(ang)
    problematic = ang < 3 * np.pi/4

    return ang, direction, problematic, negative



if __name__ == "__main__":
    obs = Obstacle(1, 2.0, 1.0, 4.1)
    # obs.swap()
    path = PathSegment(0, 0, 1, -0.1)
    path.swap()
    is_left = isLeft(obs, path)
    normal = np.array([obs.y1 - obs.y0, obs.x0 - obs.x1])/10.0
    print(is_left)
    angle, direction, problem, negative = get_angle(obs, path, is_left)
    print(problem)
    print(angle*180/np.pi)
    if direction == -1:
        plt.plot([obs.x0], [obs.y0], 'ro')
    else:
        plt.plot([obs.x1], [obs.y1], 'bo')

    if negative:
        plt.plot([obs.x0 - normal[0]], [obs.y0 - normal[1]], 'ro')
    else:
        plt.plot([obs.x0 + normal[0]], [obs.y0 + normal[1]], 'ro')
    plt.plot([obs.x0, obs.x1], [obs.y0, obs.y1], 'r')
    plt.plot([path.x0, path.x1], [path.y0, path.y1], 'b')
    plt.gca().set_aspect("equal")
    plt.show()

