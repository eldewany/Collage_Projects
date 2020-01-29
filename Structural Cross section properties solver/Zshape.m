classdef Zshape
    
    properties
        b_z(1,2)
        h_z(1,2)
        
    end
    
    methods
        
        function A=a_z(obj)
            A=(obj.b_z(1).*obj.h_z(1))+2.*(obj.b_z(2).*obj.h_z(2));
        end
        
        function Iy=iy(obj)
           
            Iy=(obj.b_z(1).*obj.h_z(1))./12.*(obj.h_z(1).^2)+2.*((obj.b_z(2).*obj.h_z(2))./12.*(obj.h_z(2).^2)+((obj.b_z(2).*obj.h_z(2)).*((obj.h_z(2)+obj.h_z(1))./2).^2));
        end
        function Iz=iz(obj)
            
            Iz=(obj.b_z(1).*obj.h_z(1))./12.*(obj.b_z(1).^2)+2.*(((obj.b_z(2).*obj.h_z(2))./12*obj.b_z(2).^2)+((obj.b_z(2).*obj.h_z(2)).*((obj.b_z(2)-obj.b_z(1))./2).^2));
        end
        function IYZ=Iyz(obj)
            IYZ=2.*(obj.b_z(2).*obj.h_z(2).*(obj.h_z(1)./2-obj.h_z(2)./2).*(obj.b_z(1)./2+obj.b_z(2)./2));
        end
    end
end

        
        
        
        
