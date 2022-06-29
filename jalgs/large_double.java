/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class large_double implements Cloneable{
	public double number;
	public int exponent;
	public boolean negative;

	public large_double(){
		number=0.0;
		exponent=0;
		negative=false;
	}

	public large_double(double number,int exponent){
		double absnumber=number;
		negative=false;
		if(absnumber<0.0){
			negative=true;
			absnumber=-number;
		}
		generate(absnumber,exponent,negative);
	}

	public large_double(double number){
		this(number,0);
	}

	public large_double(double number,int exponent,boolean negative){
		// note that the sign on number is ignored here
		generate(number,exponent,negative);
	}

	public void generate(double number,int exponent,boolean negative){
		this.negative=negative;
		double absnumber=Math.abs(number);
		if(absnumber<1.0){
			if(absnumber==0.0){
				this.number=0.0;
				this.exponent=0;
			}else{
				int added=(int)(Math.log10(absnumber))-1;
				double rem=absnumber/(Math.pow(10.0,added));
				this.number=rem;
				this.exponent=exponent+added;
			}
		}else{
			if(absnumber>=10.0){
				int added=(int)(Math.log10(absnumber));
				double rem=absnumber/(Math.pow(10.0,added));
				this.number=rem;
				this.exponent=exponent+added;
			}else{
				this.number=absnumber;
				this.exponent=exponent;
			}
		}
	}

	public void divide(double denominator){
		multiply(1.0/denominator);
	}

	public void multiply(double multiplier){
		if(multiplier==0.0){
			negative=false;
			number=0.0;
			exponent=0;
			return;
		}
		double absmultiplier=multiplier;
		if(absmultiplier<0.0){
			if(negative){
				negative=false;
			}else{
				negative=true;
			}
			absmultiplier=-multiplier;
		}
		double newnumber=absmultiplier*number;
		generate(newnumber,exponent,negative);
	}

	public void divide(large_double denominator){
		if(denominator.negative){
			if(negative){
				negative=false;
			}else{
				negative=true;
			}
		}
		double newnumber=number/denominator.number;
		exponent-=denominator.exponent;
		generate(newnumber,exponent,negative);
	}

	public void multiply(large_double multiplier){
		if(multiplier.number==0.0){
			negative=false;
			number=0.0;
			exponent=0;
			return;
		}
		if(multiplier.negative){
			if(negative){
				negative=false;
			}else{
				negative=true;
			}
		}
		double newnumber=multiplier.number*number;
		exponent+=multiplier.exponent;
		generate(newnumber,exponent,negative);
	}

	public large_double divide_copy(double divisor){
		return mult_copy(1.0/divisor);
	}

	public large_double mult_copy(double multiplier){
		if(multiplier==0.0){
			return new large_double();
		}else{
			if(multiplier<0.0){
				return new large_double(number*multiplier,exponent,!negative);
			}else{
				return new large_double(number*multiplier,exponent,negative);
			}
		}
	}

	public large_double divide_copy(large_double divisor){
		if(divisor.negative){
			return new large_double(number/divisor.number,exponent-divisor.exponent,!negative);
		}else{
			return new large_double(number/divisor.number,exponent+divisor.exponent,negative);
		}
	}

	public large_double mult_copy(large_double multiplier){
		if(multiplier.number==0.0){
			return new large_double();
		}else{
			if(multiplier.negative){
				return new large_double(number*multiplier.number,exponent+multiplier.exponent,!negative);
			}else{
				return new large_double(number*multiplier.number,exponent+multiplier.exponent,negative);
			}
		}
	}

	public void intpow(int power){
		double temp=number;
		for(int i=1;i<power;i++){
			temp*=number;
		}
		if(negative){
			if((power%2)==0){
				generate(temp,exponent*power,false);
			}else{
				generate(temp,exponent*power,true);
			}
		}else{
			generate(temp,exponent*power,false);
		}
	}

	public static large_double factorial(int val){
		large_double temp=new large_double(val);
		for(int i=(val-1);i>1;i--){
			temp.multiply(i);
		}
		return temp;
	}

	public int compare(large_double compval){
		if(number==Double.NaN || compval.number==Double.NaN) return -2;
		if(compval.number==0.0){
			if(number==0.0){
				return 0;
			}else{
				if(negative){
					return -1;
				}else{
					return 1;
				}
			}
		}
		if(number==0.0){
			if(compval.negative){
				return 1;
			}else{
				return -1;
			}
		}
		if(!negative){
			if(compval.exponent==exponent&&compval.number==number){
				if(compval.negative){
					return 1;
				}else{
					return 0;
				}
			}else{
				if(compval.exponent<exponent){
					return 1;
				}else{
					if(compval.exponent>exponent){
						return -1;
					}else{
						if(compval.number<number){
							return 1;
						}else{
							return -1;
						}
					}
				}
			}
		}else{
			if(compval.exponent==exponent&&compval.number==number){
				if(compval.negative){
					return 0;
				}else{
					return -1;
				}
			}else{
				if(compval.exponent<exponent){
					return -1;
				}else{
					if(compval.exponent>exponent){
						return 1;
					}else{
						if(compval.number<number){
							return -1;
						}else{
							return 1;
						}
					}
				}
			}
		}
	}

	public large_double clone(){
		return new large_double(number,exponent,negative);
	}

	public double get_number(){
		return number;
	}

	public int get_exponent(){
		return exponent;
	}

	public boolean is_negative(){
		return negative;
	}

	public String toString(){
		if(negative){
			return "-"+number+"E"+exponent;
		}else{
			return ""+number+"E"+exponent;
		}
	}

}
