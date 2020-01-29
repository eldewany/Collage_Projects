classdef Tshape
    
    properties
        b_t(1,2)
        h_t(1,2)
        
    end
    
    methods
        
        function A=a_t(obj)
            A=(obj.b_t(1).*obj.h_t(1))+(obj.b_t(2).*obj.h_t(2));
        end
        
        function Iy=iy(obj)
            Z_bar=((obj.b_t(1).*obj.h_t(1).*obj.h_t(1))./2+(obj.b_t(2).*obj.h_t(2)).*(obj.h_t(1)+obj.h_t(2)/2))./((obj.b_t(1).*obj.h_t(1))+(obj.b_t(2).*obj.h_t(2)));
            Iy=(obj.b_t(1).*obj.h_t(1)).*(obj.h_t(1)./2-Z_bar).^2+(obj.b_t(1)./12.*obj.h_t(1).^3)+(obj.b_t(2).*obj.h_t(2)).*(obj.h_t(1)+obj.h_t(2)./2-Z_bar).^2+(obj.b_t(2)./12.*obj.h_t(2).^3);
        
        end
        function Iz=iz(obj)
            Iz=(obj.b_t(1).*obj.h_t(1))./12.*obj.b_t(1).^2+(obj.b_t(2).*obj.h_t(2))./12*obj.b_t(2).^2;
        end
    end
end

        
        
        
        
