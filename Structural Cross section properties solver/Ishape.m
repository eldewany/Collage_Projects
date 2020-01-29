classdef Ishape
    
    properties
        b_i(1,2)
        h_i(1,2)
        
    end
    
    methods
        
        function A=a_i(obj)
            A=(obj.b_i(1).*obj.h_i(1))+2.*(obj.b_i(2).*obj.h_i(2));
        end
        
        function Iy=iy(obj)
           
            Iy=(obj.b_i(1).*obj.h_i(1))./12.*(obj.h_i(1).^2)+2.*((obj.b_i(2).*obj.h_i(2))./12.*(obj.h_i(2).^2)+(obj.b_i(2).*obj.h_i(2)).*(((obj.h_i(1)+obj.h_i(2))./2).^2));
        
        end
        function Iz=iz(obj)
            Iz=(obj.b_i(1).*obj.h_i(1))./12.*(obj.b_i(1).^2)+((obj.b_i(2).*obj.h_i(2))./12*obj.b_i(2).^2).*2;
        end
    end
end

        
        
        
        
