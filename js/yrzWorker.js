
let previousProgress = 0;

const kb_ = 8.62e-5;
const R_ = 8.314;
const alphaBCS_ = 2.14;

let alpha_ = alphaBCS_;
let Tc_ = 80.0;
let Eg_ = 2.0*alphaBCS_*kb_*Tc_;
let thc_ = Math.PI/4.0;
let Delta0_ = 0;
let x_ = 0.5;
let mu_ = -0.46;
let t0_ = 1.0;

const Nk_ = 3000;

let energy_ = [];
let dos_ = [];
let NSdos_=[];
let bin3_ = [];
let bin4_ = [];
let E1_ = -1.5*t0_;
let E2_ = 2.5*t0_;
const NE_ = 2501;
let dE_ = (E2_-E1_)/(NE_-1);

let gap_ = [];
let entropy_ = [];
let NSentropy_ = [];
let temperature_ = [];
let gamma_ = [];
let NSgamma_ = [];
let DF_ = [];
let Dt_ = [];
let lambda_ = [];
let T1_ = 2*t0_;
let T2_ = 1200*t0_;
const NT_ = 301;
let dT_ = (T2_-T1_)/(NT_-1);




let nTc_ = 80;//0.8*(NT_-1);
//let dT_ = (Tc_-T1_)/nTc_;
//let T2_ = T1_+(NT_-1)*dT_;


let t_ = [t0_,-0.3*t0_,0.2*t0_];
	
let chi_ = 0.338;
let J_ = 1.0/3.0*t_[0];




onmessage = function(event) {
    
    alpha_ = event.data.alpha;
    x_ = parseFloat(event.data.x);
    mu_ = parseFloat(event.data.mu);
    Eg_ = 1.0*3.0*(0.2 - x_);
    Delta0_ = 0.14*(1.0 - 82.6*Math.pow(x_ - 0.16, 2));
    if(Delta0_<0) Delta0_ = 0;
    Tc_ = Delta0_/alpha_/kb_;
    nTc_ = 100;
    dT_ = (Tc_-T1_)/nTc_;
    T2_ = T1_+(NT_-1)*dT_;
    /*
    this.console.log(Tc_);
    this.console.log(T2_);
    this.console.log(Delta0_);
    this.console.log(x_);
    this.console.log(mu_);
    this.console.log(Eg_);
    this.console.log(calculateDoping());
    */

    let p = calculateDoping();
    postMessage(
        {messageType: "Parameters", doping: p, Eg: Eg_}
      );

    initArrays();
    computeDos(0);
    for(let i = 0; i<NE_; i++){
        NSdos_[i] = bin3_[i]+bin4_[i];
    }
    computeDos(Delta0_);
    for(let i = 0; i<NE_; i++){
        dos_[i] = bin3_[i]+bin4_[i];
    }
    postMessage(
        {messageType: "Dos", E: energy_, N: dos_, Nn: NSdos_}
      );

    calculateDelta();  
    postMessage(
        {messageType: "Delta", D: gap_}
      );
    
    
    calculateEntropy(alpha_, Tc_);  // Send back the results.
    postMessage(
        {messageType: "Entropy", S: entropy_, Sn: NSentropy_}
      );
    calculateGamma();
    postMessage(
        {messageType: "Gamma", g: gamma_, gn: NSgamma_}
      );
    calculateDF();
    postMessage(
        {messageType: "DF", F: DF_}
      );
      
    calculateDt();
    postMessage(
        {messageType: "Dt", Dt: Dt_}
      );
      
    //calculateLvsT();
    postMessage(
        {messageType: "Lambda", L: lambda_}
      );

      
    postMessage(
        {messageType: "Done"}
      );  
    
  };
  

  function initArrays(){
    energy_ = new Array(NE_);
    for(let i = 0; i<NE_; i++){
        energy_[i] = (E1_*100+i*dE_*100)/100;
    }    
    dos_ = new Array(NE_).fill(0);
    NSdos_ = new Array(NE_).fill(0);

    bin3_ = new Array(NE_).fill(0);
    bin4_ = new Array(NE_).fill(0);

    temperature_ = new Array(NT_);
    for(let i = 0; i<NT_; i++){
        temperature_[i] = T1_+i*dT_;
    }
    gap_ = new Array(NT_).fill(0);
    entropy_ = new Array(NT_).fill(0);
    NSentropy_ = new Array(NT_).fill(0);
    //gamma_ = new Array(NT_).fill(0);
    //NSgamma_ = new Array(NT_).fill(0);
    DF_ = new Array(NT_).fill(0);
    Dt_ = new Array(NT_).fill(0);
    lambda_ = new Array(NT_).fill(0);
}


function computeDos(Delta){

	//EvHs_ = 0; //xi0(PI_,0,x)-4.0*tp(x)*cos(PI_)*cos(0)-2.0*tpp(x)*(cos(2.0*PI_)+cos(2.0*0));
    for(let i = 0; i<NE_; i++){
        bin3_[i] = 0;
        bin4_[i] = 0;
    }
    
    let Ep = 0;
    let Em = 0;
    let Eps = 0;
    let Ems = 0;
    let up2 = 1.0;
	let	um2 = 1.0;
	let	vp2 = 0;
	let	vm2 = 0;
	
    let kx = 0;
    let ky = 0;
	let k1 = 0;
	let k2 = Math.PI;
	let dk = (k2-k1)/(Nk_-1);
	//let V,Vp,Vm;
    let gtWp = 0;
    let gtWm = 0;
	let gtVar = 0;
    let indexp = 0;
    let indexm = 0;
    let indexp2 = 0;
    let indexm2 = 0;
	let mut = mu_*t0_;
    let Egt = Eg_*t0_;
    let Deltat = Delta*t0_;
    
	
	let progress = 1;
	
	for(let i = 0; i<Nk_; i++){
		kx = k1+i*dk;
		for(let j = 0; j<Nk_; j++){
            ky = k1+j*dk;
            
            Dk = pseudogap(kx,ky,Deltat);
            Ep = Ek(kx,ky,x_,mut,Egt,1);
			Em = Ek(kx,ky,x_,mut,Egt,-1);
			gtVar = gt(x_);
			gtWp = gtVar*Wk(kx,ky,x_,mut,Egt,1);
			gtWm = gtVar*Wk(kx,ky,x_,mut,Egt,-1);
            
            
            if(Delta==0){
				Eps = Ep;
				Ems = Em;
				up2 = 1.0;
				um2 = 1.0;
				vp2 = 0;
				vm2 = 0;
			}else{	
				Eps = Math.sqrt(Math.pow(Ep,2)+Math.pow(Dk,2));
				Ems = Math.sqrt(Math.pow(Em,2)+Math.pow(Dk,2));
				up2 = 0.5*(1.0+Ep/Eps);
				um2 = 0.5*(1.0+Em/Ems);
				vp2 = 0.5*(1.0-Ep/Eps);
				vm2 = 0.5*(1.0-Em/Ems);
			}

            indexp = getIndex(Eps);
			indexp2 = getIndex(-Eps);
			indexm = getIndex(Ems);
            indexm2 = getIndex(-Ems);
            
            tally(indexp,gtWp*up2,bin3_);
			tally(indexm,gtWm*um2,bin4_);
			tally(indexp2,gtWp*vp2,bin3_);
			tally(indexm2,gtWm*vm2,bin4_);

			//W = 1.0;
			
			//V = dxi(kx,ky,x);
			//V = pow((cos(kx)-cos(ky))/2.0,2); //c-axis
			
			//indexp = getIndex(Ep);
			//indexm = getIndex(Em);
			
						
			//yrz
			//tally(indexp,gtWp,bin3_);
			//tally(indexm,gtWm,bin4_);
			
			
			/*
			tally(indexp,pow(V*gtWp,2),bin_);
			tally(indexm,pow(V*gtWm,2),bin2_);
			if(indexp==indexm){
				tally(indexp,2.0*V*V*gtWp*gtWm,bin_);
			}	
			
			V = pow(dxi(kx,ky,x),2)*d2xidy2(kx,ky,x)-dxi(kx,ky,x)*dxidy(kx,ky,x)*d2xidydx(kx,ky,x); //hall effect
			tally(indexp,V*pow(gtWp,3),bin5_);
			tally(indexm,V*pow(gtWm,3),bin6_);
			if(indexp==indexm){
				tally(indexp,3.0*V*(pow(gtWp,2)*gtWm+gtWp*pow(gtWm,2)),bin5_);
            }
            */
			
			
        }
        /*
		if(Math.floor(100*i/(progress*Nk_))==1){
            progress++;
            postMessage(
                {messageType: "Progress", data: progress}
                );
        }
        */
	}
	for(let i = 0; i<NE_; i++){
		//bin_[i] = bin_[i]/Nk_/Nk_/dE_;
		//bin2_[i] = bin2_[i]/Nk_/Nk_/dE_;
		bin3_[i] = bin3_[i]/Nk_/Nk_/dE_;
        bin4_[i] = bin4_[i]/Nk_/Nk_/dE_;
        
        //console.log(dos_[i]);
 		//bin5_[i] = bin5_[i]/Nk_/Nk_/dE_;
		//bin6_[i] = bin6_[i]/Nk_/Nk_/dE_;
	}
	
}


function tally(index, amount, bin){
	if(index<NE_ && !(index<0)){
		bin[index] = bin[index]+amount;
	}
}


function getIndex(energy){
	
	let d = (energy - E1_)/dE_;
	
	let index = Math.floor(d);
	
	if((d-index)>=0.5){	
		index++;
    }
    
	return index;
	
}


function gt(x){
	//return 2.0*0.2/(1.0+0.2);
	return 2.0*x/(1.0+x);
}

function gs(x){
	//return 4.0/pow(1.0+0.2,2);
	return 4.0/Math.pow(1.0+x,2);
}

function t(x){
	return gt(x)*t_[0]+3.0/8.0*gs(x)*J_*chi_;
}	

function tp(x){
	return gt(x)*t_[1];
}

function tpp(x){
	return gt(x)*t_[2];
}

function tperp(kx, ky){
	return 0.0;//*t_[0]*pow((cos(kx)-cos(ky))/2.0,2);
}

function xi0(kx, ky, x){
	return -2.0*t(x)*(Math.cos(kx)+Math.cos(ky));
}

function xi(kx, ky, mu, x){
    //yrz with mu relative to vHs.
	let energy = xi0(kx,ky,x)-4.0*tp(x)*Math.cos(kx)*Math.cos(ky)-2.0*tpp(x)*(Math.cos(2.0*kx)+Math.cos(2.0*ky))-mu-tperp(kx,ky);
	return energy;
	
}

function pseudogap(kx, ky, Eg){
	return Eg/2.0*(Math.cos(kx)-Math.cos(ky));
}

//Renormalized dispersion.
//sign must be -1 or +1
function Ek(kx, ky, x, mu, Eg, sign){
	let xiVar = xi(kx,ky,mu,x);
	let xi0Var = xi0(kx,ky,x);
	return 0.5*(xiVar-xi0Var)+sign*Math.sqrt(Math.pow(0.5*(xiVar+xi0Var),2)+Math.pow(pseudogap(kx,ky,Eg),2));
}

//Weight.
//sign must be -1 or +1
function Wk(kx, ky, x, mu, Eg, sign){
	let xiVar = xi(kx,ky,mu,x);
	let xi0Var = xi0(kx,ky,x);
	return 0.5*(1.0+sign*(xiVar+xi0Var)/2.0/Math.sqrt(Math.pow((xiVar+xi0Var)/2.0,2)+Math.pow(pseudogap(kx,ky,Eg),2)));
}





  function calculateGamma(){
    for(let i = 0; i<NT_-1; i++){
        gamma_[i]=(entropy_[i+1]-entropy_[i])/dT_;
        NSgamma_[i]=(NSentropy_[i+1]-NSentropy_[i])/dT_;
    }
}

function calculateDt(){
    let x = 0;
    for(let i = 0; i<NT_; i++){
        x = 0;
        if(temperature_[i]<Tc_)
        x = Math.sqrt(DF_[i]/DF_[0])-(1.0-Math.pow(temperature_[i]/Tc_,2));
        Dt_[i] = x;
    }
}

function calculateDF(){
    let sum = 0;
    for(let i = nTc_-1; i>=0; i--){
        sum+=0.5*(entropy_[i]-NSentropy_[i]+entropy_[i+1]-NSentropy_[i+1]);
        DF_[i]=-sum*dT_;
    }
}



  function calculateEntropy(alpha, Tc) {
   
    // Calculate the progress percentage.
    let progress = 0;//Math.round(i/list.length*100);  // Only send a progress update if the progress has changed
    // at least 1%.
    let zeroGapReached = false;
    computeDos(0);
    for(let i = 0; i<NT_; i++){
        if(temperature_[i]==0){
            NSentropy_[i] = 0;
        }else{
            NSentropy_[i] = S(temperature_[i]);
        }    

    }


    for(let i = 0; i<NT_; i++){
        progress = Math.floor(i/(NT_-1)*100);
        if(!zeroGapReached){ 
            computeDos(gap_[i]/1000);
            if(gap_[i]==0){
                zeroGapReached = true;
            }
        }
        if(temperature_[i]==0){
            entropy_[i] = 0; 
        }else{
            entropy_[i] = S(temperature_[i]);
        }
        lambda_[i] = calculateSFD(gap_[i]/1000,temperature_[i]);    

        if (progress != previousProgress) {
            postMessage(
            {messageType: "Progress", data: progress}
            );
            previousProgress = progress;
        }
    }   
  }



 function S(T){
	let mu = 0;
    let sum = 0;
    let fw = 0;
    let f = fermi(E1_,mu,T);
    fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f)
    if(!isNaN(fw)){
		sum = 0.5*fw*(bin3_[0]+bin4_[0]);
	}
    f = fermi(E1_+(NE_-1)*dE_,mu,T);
    fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f)
    if(!isNaN(fw)){
		sum += 0.5*fw*(bin3_[NE_-1]+bin4_[NE_-1]);
	}
	for(let i = 1; i<NE_-1; i++){
        f = fermi(E1_+i*dE_,mu,T);
        fw = f*Math.log(f)+(1.0-f)*Math.log(1.0-f)
        if(!isNaN(fw)){
			sum += fw*(bin3_[i]+bin4_[i]);
		}
	}
	sum = -2.0*1000*R_*sum*dE_; //mJ/mol/K^2
	return sum;
}



  





function calculateDelta(){
    for(let i=0; i<NT_; i++){
        gap_[i] = 1E3*delta(alpha_,Tc_,temperature_[i]);
    }
} 


function delta(alpha,Tc,T){
    let a = 7.0/5.0;
    let DCoverC = 0.955;
    let y = a*DCoverC*(Tc/T-1.0);
    if(y<0) return 0;
    return alpha*kb_*Tc*Math.tanh(Math.PI/alphaBCS_*Math.sqrt(y));
}


function fermi(E, mu, T){
    if(T==0){
        if(E<mu){
            return 1.0;
        }
        else if(E==mu){ 
            return 0.5;
        }
        else if(E>mu){
            return 0;
        }
    }
    return 1.0/(Math.exp((E-mu)/kb_/T)+1.0);
}


function dfdt(t,Delta,T){
    let x = Math.sqrt(Math.pow(t,2)+Math.pow(Delta,2));
    let y = Math.exp(x/kb_/T);
    return -t*y/(Math.pow(y+1,2)*x*kb_*T);
}

function dfdE(E,T){
    let y = Math.exp(E/kb_/T);
    let y2 = -y/(Math.pow(y+1,2)*kb_*T);
    if(isNaN(y2)) return 0;
    return y2;
}

/*
function calculateLvsT(){
    let th1 = 0;
    let th2 = Math.PI/4.0;
    const Nth = 200;
    let dth = (th2-th1)/(Nth-1);
    let th = 0;
    let d0 = 0;
    let pg = 0;
    for(let i = 0; i<NT_; i++){
        d0 = delta(alpha_,Tc_,temperature_[i]);
        for(let j = 0; j<Nth; j++){
            th = th1+j*dth;
            //pg = Eg_*(1.0-th/thc_);
            pg = 0;
            if(th<thc_)
                pg = Eg_*Math.cos(2.0*Math.PI*th/4.0/thc_);
           
            lambda_[i] += calculateLambda(0.02,d0*Math.cos(2.0*th),pg, temperature_[i]);
        }
        lambda_[i] = lambda_[i]*dth/Math.PI*4.0;  
        //console.log(temperature_[i]);  
    }
}






function calculateLambda(mu, Delta, pg, T){
    E1_ = -0.1;
    E2_ = -E1_;
    let NEL = 2*NE_;
    dE_ = (E2_-E1_)/(NEL-1);
    let E = 0;
    let a = 0;
    let b = 0;
    let c = 0;
    let d = 0;
    let x = 0;
    let sum = 0;

    E = E1_;
    if(Math.abs(E)>=pg){
        if(E==0){
            sum+=-0.5*2.0*mu*(Math.exp(Delta/kb_/T)/Math.pow(Math.exp(Delta/kb_/T)+1.0,2)/kb_/T)+0.5*mu/Delta*(1.0-2.0*fermi(Delta,0,T));
        }else{
        x = Math.sqrt(Math.pow(E,2)+Math.pow(Delta,2));
        a = Math.pow(Delta,2)/x;
        b = mu*Math.pow(Delta,2)/E/x;
        c = 2.0*dfdt(E,Delta,T);
        if(isNaN(c)) c=0;
        d = E*(1.0-2.0*fermi(x,0,T))/Math.pow(x,2);
        sum+=0.5*(a+b)*(c+d);
        }
    }

        E = E2_;
        if(Math.abs(E)>=pg){
        if(E==0){
            sum+=-0.5*2.0*mu*(Math.exp(Delta/kb_/T)/Math.pow(Math.exp(Delta/kb_/T)+1.0,2)/kb_/T)+0.5*mu/Delta*(1.0-2.0*fermi(Delta,0,T));
        }else{
        x = Math.sqrt(Math.pow(E,2)+Math.pow(Delta,2));
        a = Math.pow(Delta,2)/x;
        b = mu*Math.pow(Delta,2)/E/x;
        c = 2.0*dfdt(E,Delta,T);
        if(isNaN(c)) c=0;
        d = E*(1.0-2.0*fermi(x,0,T))/Math.pow(x,2);
        sum+=0.5*(a+b)*(c+d);
        }    
        }


    for(let i = 1; i<NEL-1; i++){
        E = E1_+i*dE_;
        if(Math.abs(E)>=pg){
        if(E==0){
            sum+=-2.0*mu*(Math.exp(Delta/kb_/T)/Math.pow(Math.exp(Delta/kb_/T)+1.0,2)/kb_/T)+mu/Delta*(1.0-2.0*fermi(Delta,0,T));
        }else{
        x = Math.sqrt(Math.pow(E,2)+Math.pow(Delta,2));
        a = Math.pow(Delta,2)/x;
        b = mu*Math.pow(Delta,2)/E/x;
        c = 2.0*dfdt(E,Delta,T);
        if(isNaN(c)) c=0;
        d = E*(1.0-2.0*fermi(x,0,T))/Math.pow(x,2);
        if(isNaN(d)) console.log(d);
        sum+=(a+b)*(c+d);
        }
        }
    }
    return sum*=dE_;
}
*/


function calculateSFD(delta, T) {

	if (delta == 0) {
		return 0;
	}


    let E = 0;
    let Ep = 0;
    let Em = 0;
    let Eps = 0;
    let Ems = 0;
    let up = 0;
    let um = 0;
    let vp = 0;
    let vm = 0;
    let Dk = 0;
    let dfEps = 0;
    let dfEms = 0;
    let fEps = 0;
    let fEms = 0;
    let A = 0;
    let B = 0;
    let C = 0;

    let kx = 0;
    let ky = 0;
	let k1 = 0;
	let k2 = Math.PI;
	let dk = (k2 - k1) / (Nk_ - 1);
	let V = 0;
    let Wp = 0; 
    let Wm = 0;

	let vxp = 0;
	let vxm = 0;
	let vyp = 0;
	let vym = 0;
	let vyyp = 0;
	let vyym = 0;
	let vyxp = 0;
    let vyxm = 0;
    

    let xix = 0;
    let xiy = 0;
    let xiyy = 0;
    let xiyx = 0;

    let xi0x = 0;
    let xi0y = 0;
    let xi0yy = 0;
    let xi0yx = 0;
    
    let pgx = 0;
    let pgy = 0;
    let pgyy = 0;
    let pgyx = 0;

    let a1 = 0;
    let a2 = 0;
    let a3 = 0;
    let a4 = 0;
    let a5 = 0;
    let a6 = 0;
    let a7 = 0;
    let a8 = 0;
    let a9 = 0;
    let a10 = 0;


	let Bp = 0;
	let Bm = 0;
	let sigmap = 0;
	let sigmam = 0;

	let w = 0;

	let Deltat = delta*t0_;

	let mut = mu_*t0_;
	let Egt = Eg_*t0_;

	let progress = 1;
    
    let xiVar = 0;
    let xi0Var = 0;
    let pg = 0;
    let y = 0;

	let rhos = 0;

	for (let i = 0; i<Nk_; i++) {
		kx = k1 + i*dk;
		for (let j = 0; j<Nk_; j++) {
			ky = k1 + j*dk;


			Dk = pseudogap(kx, ky, Deltat);
			//getBandValues(kx, ky, x, mut, Egt, Ep, Em, Wp, Wm);
			//getDerivs(kx, ky, x, mut, Egt, vxp, vxm, vyp, vym, vyyp, vyym, vyxp, vyxm);

            xiVar = xi(kx, ky, mut, x_);
	        xi0Var = xi0(kx, ky, x_);
	        pg = pseudogap(kx, ky, Egt);
	        y = Math.sqrt(Math.pow(0.5*(xiVar + xi0Var), 2) + Math.pow(pg, 2));
	        Ep = 0.5*(xiVar - xi0Var) + y;
	        Em = 0.5*(xiVar - xi0Var) - y;
	        Wp = 0.5*(1.0 + (xiVar + xi0Var) / 2.0 / y);
	        Wm = 0.5*(1.0 - (xiVar + xi0Var) / 2.0 / y);

            xix = dxi(kx, ky, x_);
            xiy = dxidy(kx, ky, x_);
            xiyy = d2xidy2(kx, ky, x_);
            xiyx = d2xidydx(kx, ky, x_);

            xi0x = dxi0dx(kx, ky, x_);
            xi0y = dxi0dy(kx, ky, x_);
            xi0yy = d2xi0dy2(kx, ky, x_);
            xi0yx = d2xi0dxdy(kx, ky, x_);
            
            pgx = dPGdx(kx, ky, Egt);
            pgy = dPGdy(kx, ky, Egt);
            pgyy = d2PGdy2(kx, ky, Egt);
            pgyx = d2PGdxdy(kx, ky, Egt);

            a1 = 0.5*(xiy - xi0y);
            a2 = 0.5*(xiVar + xi0Var)*(xiy + xi0y) + 2.0*pg*pgy;
            vyp = a1 + 0.5*a2 / y;
            vym = a1 - 0.5*a2 / y;

            a3 = 0.5*(xix - xi0x);
            a4 = 0.5*(xiVar + xi0Var)*(xix + xi0x) + 2.0*pg*pgx;
            vxp = a3 + 0.5*a4 / y;
            vxm = a3 - 0.5*a4 / y;

            a5 = 0.5*(xiyy - xi0yy);
            a6 = 0.25*Math.pow(xiy + xi0y, 2) + 0.25*(xiVar + xi0Var)*(xiyy + xi0yy) + Math.pow(pgy, 2) + pg*pgyy;
            vyyp = a5 - 0.25*Math.pow(a2, 2) / Math.pow(y, 3) + a6 / y;
            vyym = a5 + 0.25*Math.pow(a2, 2) / Math.pow(y, 3) - a6 / y;

            a7 = 0.5*(xiyx - xi0yx);
            a8 = 0.5*(xiVar + xi0Var)*(xix + xi0x) + 2.0*pg*pgx;
            a9 = -0.25*a2*a8 / Math.pow(y, 3);
            a10 = 0.25*(xix + xi0x)*(xiy + xi0y) + 0.25*(xiVar + xi0Var)*(xiyx + xi0yx) + pgx*pgy + pg*pgyx;
            vyxp = a7 +a9 + a10 / y;
            vyxm = a7 -a9 - a10 / y;

			if (delta == 0) {
				Eps = Ep;
				Ems = Em;
				//up = 1.0;
				//um = 1.0;
				//vp = 0;
				//vm = 0;
			}
			else {
				Eps = Math.sqrt(Math.pow(Ep, 2) + Math.pow(Dk, 2));
				Ems = Math.sqrt(Math.pow(Em, 2) + Math.pow(Dk, 2));
				up = Math.sqrt(0.5*(1.0 + Ep / Eps));
				um = Math.sqrt(0.5*(1.0 + Em / Ems));
				vp = Math.sqrt(0.5*(1.0 - Ep / Eps));
				vm = Math.sqrt(0.5*(1.0 - Em / Ems));
			}

			V = dxi(kx, ky, x_);
			dfEps = dfdE(Eps, T);
			dfEms = dfdE(Ems, T);
			fEps = fermi(Eps, 0, T);
			fEms = fermi(Ems, 0, T);

			A = Wp*Wp*Math.pow(up*vp, 2)*(2.0*dfEps + (1.0 - 2.0*fEps) / Eps);
			B = Wm*Wm*Math.pow(um*vm, 2)*(2.0*dfEms + (1.0 - 2.0*fEms) / Ems);
			//C = 4.0*Wp*Wm*up*vp*um*vm*((fEps - fEms) / (Eps - Ems) + (1.0 - fEps - fEms) / (Eps + Ems));
			//rhos += V*V*(A + B + C);
			
			rhos += vxp*vxp*A + vxm*vxm*B;

			

		}
		if (100 * i / (progress*Nk_) == 1) {
			progress++;
		}
	}

	
	//return 16.0*PI_ / t0_*gt(x)*gt(x)*rhos / Nk_ / Nk_;
	return 16.0*Math.PI / t0_*rhos / Nk_ / Nk_;

}



function dxi(kx, ky, x){
	return 2.0*t(x)*Math.sin(kx)+4.0*tp(x)*Math.sin(kx)*Math.cos(ky)+4.0*tpp(x)*Math.sin(2.0*kx); //yrz
}

function dxidy(kx, ky, x){
	return 2.0*t(x)*Math.sin(ky)+4.0*tp(x)*Math.cos(kx)*Math.sin(ky)+4.0*tpp(x)*Math.sin(2.0*ky); //yrz
}

function d2xidy2(kx, ky, x){
	return 2.0*t(x)*Math.cos(ky)+4.0*tp(x)*Math.cos(kx)*Math.cos(ky)+8.0*tpp(x)*Math.cos(2.0*ky); //yrz
}

function d2xidydx(kx, ky, x){
	return -4.0*tp(x)*Math.sin(kx)*Math.sin(ky); //yrz
}

function dxi0dx(kx, ky, x) {
	return 2.0*t(x)*Math.sin(kx);
}

function dxi0dy(kx, ky, x) {
	return 2.0*t(x)*Math.sin(ky);
}

function d2xi0dy2(kx, ky, x) {
	return 2.0*t(x)*Math.cos(ky);
}

function d2xi0dxdy(kx, ky, x) {
	return 0;
}

function dPGdx(kx, ky, Eg) {
	return -Eg / 2.0 * Math.sin(kx);
}

function dPGdy(kx, ky, Eg) {
	return Eg / 2.0*Math.sin(ky);
}

function d2PGdy2(kx, ky, Eg) {
	return Eg / 2.0*Math.cos(ky);
}

function d2PGdxdy(kx, ky, Eg) {
	return 0;
}



function calculateDoping() {

	let kx, ky;
	let k1 = 0;
	let k2 = Math.PI;
	let tol = 0.01;
	let N = 3000;
	let dk = (k2 - k1) / (N - 1);
	let Ep, Em;
	let gtVar, gtWp, gtWm;

	let holes = 0;
	let electrons = 0;
	let area;
	let p;


	for (let i = 0; i<N; i++) {
		kx = k1 + i*dk;
		for (let j = 0; j<N; j++) {
			ky = k1 + j*dk;
			Ep = Ek(kx, ky, x_, mu_*t0_, Eg_*t0_, 1.0);
			Em = Ek(kx, ky, x_, mu_*t0_, Eg_*t0_, -1.0);

			if (Eg_ == 0) {
				Em = xi(kx, ky, mu_*t0_, x_);
				Ep = 0;
			}
			if (Em>0) {
				holes += 1;
			}
			if (Ep<0) {
				electrons += 1;
			}
		}

	}

	area = (holes - electrons)*Math.pow(dk, 2);
	p = 2.0*area / Math.pow(Math.PI, 2);
	if (Eg_ == 0) {
		p -= 1.0;
	}
	else {
		//p /= 2.0; //for AFCore only
	}
	
	return p;

}
