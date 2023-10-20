import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List


@dataclass
class Obstacle:
    x0: float
    y0: float
    x1: float
    y1: float

    def draw(self, ax):
        ax.plot([self.x0, self.x1], [self.y0, self.y1], 'k', linewidth=2)


@dataclass
class Path:
    x: np.ndarray
    y: np.ndarray
    angle: np.ndarray = field(default_factory=lambda: np.zeros((0,)))

    def draw(self, ax):
        ax.plot(self.x, self.y, 'b', linewidth=2)
        for i in range(self.angle.shape[0]):
            ax.arrow(self.x[i], self.y[i], -0.5 * np.sin(self.angle[i]),
                     0.5 * np.cos(self.angle[i]), width=0.02, color='r')

    def __post_init__(self):
        self.angle = np.zeros((self.x.shape[0],))
        self.angle[0] = np.arctan2(self.y[1] - self.y[0], self.x[1] - self.x[0])
        self.angle[-1] = np.arctan2(self.y[-1] - self.y[-2],
                                    self.x[-1] - self.x[-2])
        for i in range(1, self.x.shape[0] - 1):
            self.angle[i] = np.arctan2(self.y[i + 1] - self.y[i - 1],
                                       self.x[i + 1] - self.x[i - 1])


@dataclass
class Corridor:
    path: Path
    obstacles: List[Obstacle]
    x: np.ndarray = field(default_factory=lambda: np.zeros((0,)))
    y: np.ndarray = field(default_factory=lambda: np.zeros((0,)))

    def draw(self, ax):
        ax.plot(self.x, self.y, 'b', linewidth=2)

    def __post_init__(self):
        self.x = np.zeros((self.path.x.shape[0],))
        self.y = np.zeros((self.path.x.shape[0],))
        us = np.ones((4, len(obstacles)))
        us[0, :] = np.inf
        us[1, :] = -np.inf
        for i in range(self.path.x.shape[0]):
            x0 = self.path.x[i]
            y0 = self.path.y[i]
            x1 = self.path.x[i] - np.sin(self.path.angle[i])
            y1 = self.path.y[i] + np.cos(self.path.angle[i])
            min_dist = np.inf
            obs_min = None
            u_o = None
            for obs_id, obs in enumerate(self.obstacles):
                x, y, t, u = self._intersect(x0, y0, x1, y1, obs)
                if x is not None:
                    if t < min_dist:
                        obs_min = obs_id
                        min_dist = t
                        u_o = u
                        self.x[i] = x
                        self.y[i] = y
            if u_o > us[1, obs_min]:
                us[1, obs_min] = u_o
                us[3, obs_min] = i
            if u_o < us[0, obs_min]:
                us[0, obs_min] = u_o
                us[2, obs_min] = i

        for i in range(len(obstacles)):
            self.x[int(us[2, i])] = self.obstacles[i].x0
            self.y[int(us[2, i])] = self.obstacles[i].y0
            self.x[int(us[3, i])] = self.obstacles[i].x1
            self.y[int(us[3, i])] = self.obstacles[i].y1


        # for i in range(len(self.obstacles * 2)):
        #     index = self.path.x.shape[0] + i
        #     if i % 2 == 0:
        #         self.x[index] = self.obstacles[i // 2].x0
        #         self.y[index] = self.obstacles[i // 2].y0
        #     else:
        #         self.x[index] = self.obstacles[i // 2].x1
        #         self.y[index] = self.obstacles[i // 2].y1

    def _intersect(self, x0, y0, x1, y1, obs):
        x2 = obs.x0
        y2 = obs.y0
        x3 = obs.x1
        y3 = obs.y1
        denom = (x0 - x1) * (y2 - y3) - (y0 - y1) * (x2 - x3)
        if abs(denom) < 1e-6:
            return None, None, None, None
        t = ((x0 - x2) * (y2 - y3) - (y0 - y2) * (x2 - x3)) / denom
        u = -((x0 - x1) * (y0 - y2) - (y0 - y1) * (x0 - x2)) / denom
        # if t > 0 and t < 1:
        if 0 < u < 1:
            return x0 + t * (x1 - x0), y0 + t * (y1 - y0), t, u
        else:
            return None, None, None, None

    def _intersect_full(self, x0, y0, x1, y1, obs):
        x2 = obs.x0
        y2 = obs.y0
        x3 = obs.x1
        y3 = obs.y1
        denom = (x0 - x1) * (y2 - y3) - (y0 - y1) * (x2 - x3)
        if abs(denom) < 1e-6:
            return None, None
        t = ((x0 - x2) * (y2 - y3) - (y0 - y2) * (x2 - x3)) / denom
        u = -((x0 - x1) * (y0 - y2) - (y0 - y1) * (x0 - x2)) / denom
        if 0 < u < 1 and 0 < t < 1:
            return x0 + t * (x1 - x0), y0 + t * (y1 - y0), t
        else:
            return None, None


if __name__ == "__main__":
    fig, ax = plt.subplots()
    x = np.linspace(0, 10, 20)
    y = np.exp(-x)*2 + np.sin(x) * 0.1 + x / 10
    path = Path(x, y)
    obstacle = Obstacle(2, 0.8, 4, 0.8)
    obstacle3 = Obstacle(2, 4.0, 2, 0.8)
    obstacle2 = Obstacle(-2, 4.0, 10, 4.0)
    obstacles = [obstacle, obstacle2, obstacle3]
    for ob in obstacles:
        ob.draw(ax)
    corridor = Corridor(path, obstacles)
    path.draw(ax)
    corridor.draw(ax)
    ax.set_aspect('equal')
    plt.show()
