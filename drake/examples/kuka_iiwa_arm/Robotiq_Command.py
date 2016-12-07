import os
import sys
import lcm

sys.path.append(os.path.abspath("../robotiq-driver"))
import CModelSimpleController
from robo_lcm import lcmt_robotiq_output as outputMsg

class Robotiq_Command(object):
    def __init__(self):
        self.command = outputMsg()
        self.lc = lcm.LCM()
        self.command_topic = outputMsg()
        self.command = outputMsg()

    def gripperConnect(self):
      CModelSimpleController.publisher()
      pass
    def gripperOpen(self):
      self.command = CModelSimpleController.genCommand('o',self.command)
      self.lc.publish("command_topic",self.command.encode())
    def gripperClose(self):
      self.command = CModelSimpleController.genCommand('c',self.command)
      self.lc.publish("command_topic",self.command.encode())
    def gripperReset(self):
      self.command = CModelSimpleController.genCommand('r',self.command)
      self.lc.publish("command_topic",self.command.encode())
    def gripperActivate(self):
      self.command = CModelSimpleController.genCommand('a',self.command)
      self.lc.publish("command_topic",self.command.encode())