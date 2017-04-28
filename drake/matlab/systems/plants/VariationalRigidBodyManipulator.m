classdef VariationalRigidBodyManipulator < DrakeSystem
    %This class implements the Scalable Variational Integrator algorithm
    
    properties
        manip
        timestep
        twoD = false
        dirty = true
    end
    
    methods
        function obj = VariationalRigidBodyManipulator(manipulator_or_urdf_filename,timestep,options)
            if (nargin<3)
                options=struct();
            end
            if ~isfield(options,'twoD')
                options.twoD = false;
            end
            
            typecheck(timestep,'double');
            sizecheck(timestep,1);
            
            if isempty(manipulator_or_urdf_filename) || ischar(manipulator_or_urdf_filename)
                % then make the corresponding manipulator
                w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
                if options.twoD
                    manip = PlanarRigidBodyManipulator(manipulator_or_urdf_filename,options);
                else
                    manip = RigidBodyManipulator(manipulator_or_urdf_filename,options);
                end
                warning(w);
            else
                manip = manipulator_or_urdf_filename;
            end
            
            typecheck(manip,'RigidBodyManipulator');
            obj = obj@DrakeSystem(0,manip.getNumStates(),manip.getNumInputs(),manip.getNumOutputs(),manip.isDirectFeedthrough(),manip.isTI());
            obj.manip = manip;
            if isa(manip,'PlanarRigidBodyManipulator')
                obj.twoD = true;
            end
            
            obj.timestep = timestep;
            obj = compile(obj);
        end
        
        function model = compile(model)
            model.manip = model.manip.compile();
            
            model = setNumDiscStates(model,model.manip.getNumContStates());
            model = setNumInputs(model,model.manip.getNumInputs());
            
            model = setInputLimits(model,model.manip.umin,model.manip.umax);
            model = setInputFrame(model,getInputFrame(model.manip));
            model = setStateFrame(model,getStateFrame(model.manip));

            model.dirty = false;
        end
    
        function [xdn,df] = update(obj,t,x,u)
            h = obj.timestep;
            Nq = length(x)/2;
            Nv = Nq;
            
            q0 = x(1:Nq);
            v0 = x(Nq + (1:Nv));
            
            M = manipulatorDynamics(obj.manip, q0, zeros(Nv,1));
            p0 = M*v0;
            
            q1 = q0 + h*v0;
            r = ones(size(q0));
            while max(abs(r)) > 1e-6
                r = MidpointDEL(obj,p0,q0,q1);
                dq = h*(M\r);
                alpha = 1;
                r2 = r'*r;
                rnew2 = r2+1;
                while rnew2 > r2
                    q1new = q1 + alpha*dq;
                    rnew = MidpointDEL(obj,p0,q0,q1new);
                    rnew2 = rnew'*rnew;
                    alpha = alpha/2;
                end
                q1 = q1new;
            end
            
            [D1L,D2L] = obj.LagrangianDerivs((q0+q1)/2,(q1-q0)/h);
            M = manipulatorDynamics(obj.manip, q1, zeros(Nv,1));
            p1 = -(h/2)*D1L + D2L;
            v1 = M\p1;
            
            xdn = [q1; v1];
        end
        
        function r = MidpointDEL(obj,p0,q0,q1)
            h = obj.timestep;
            [D1L,D2L] = obj.LagrangianDerivs((q0+q1)/2,(q1-q0)/h);
            %r = (h/2)*D1L1 + D2L1 + (h/2)*D1L2 - D2L2;
            r = p0 + (h/2)*D1L - D2L;
        end
        
        function [D1L,D2L] = LagrangianDerivs(obj,q,v)
            Nq = length(q);
            Nv = length(v);
            [M,G,~,dM] = manipulatorDynamics(obj.manip, q, zeros(Nv,1));
            dM = reshape(dM,Nq*Nq,Nq+Nv);
            dMdq = dM(:,1:Nq);
            D1L = 0.5*dMdq'*kron(v,v) - G;
            D2L = M*v;
        end
        
    end
    
end

