"""Follow-up derivative-bound alternatives; base depth implementation unchanged."""
from flint import arb
from arb_deep_network import DeepPoisson,activate,contract
from arb_certification import finite


class RefinedDeepPoisson(DeepPoisson):
    def __init__(self,fixture,label,sha,strategy='value_centered',bits=103):
        if strategy not in ('value_centered','centered_lite','split'):raise ValueError('Unknown refinement')
        super().__init__(fixture,label,sha,'chain',bits)
        self.strategy=strategy;self.center_passes=0;self.bound_subcells=0

    def network_jet(self,x,order):
        self.guard();finite(x,'refined input')
        if type(order) is not int or not 0<=order<=5:raise ValueError('Unsupported derivative order')
        if self.strategy=='split' or x.rad().upper()<=arb(1)/10**20:
            return super().network_jet(x,order)
        self.forward_calls+=1;self.center_passes+=1
        work=max(1,order);center_order=0 if self.strategy=='value_centered' else work
        state=[[x,arb(1)]+[arb(0)]*(work-1)]
        centers=[[x.mid()]+([arb(1)]+[arb(0)]*(center_order-1) if center_order else [])]
        radius=x.rad().upper()
        def tighten(values,mid):
            if self.strategy=='centered_lite':return contract(values,mid,radius)
            result=list(values);error=finite(radius*result[1].abs_upper(),'value-centered radius').upper()
            result[0]=finite(result[0].intersection(mid[0]+arb(0,error)),'value-centered intersection')
            return result
        for layer in range(len(self.matrices)):
            affine=self._affine(layer,state,work);ca=self._affine(layer,centers,center_order)
            affine=[tighten(z,c) for z,c in zip(affine,ca)]
            if layer<len(self.matrices)-1:
                state=[activate(z,'chain') for z in affine]
                centers=[activate(z,'chain') for z in ca]
                state=[tighten(z,c) for z,c in zip(state,centers)]
            else:state=affine
        return state[0][:order+1]

    def _split_bound(self,x,method):
        self.guard();finite(x,'derivative subcell input')
        if self.strategy!='split' or x.rad().upper()<=arb(1)/10**20:return method(x)
        midpoint=x.mid();left=x.lower().union(midpoint);right=midpoint.union(x.upper())
        self.bound_subcells+=2
        return finite(method(left).union(method(right)),'subdivided derivative union')

    def r_second(self,x):return self._split_bound(x,super().r_second)
    def r_third(self,x):return self._split_bound(x,super().r_third)
