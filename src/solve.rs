use serde::{Serialize, Deserialize};

use axgeom::*;
use axgeom::Axis;
use std::collections::BTreeMap;
use axgeom::ordered_float::*;
use duckduckgeo::grid;

    
use core::cmp::Ordering;
#[derive(Eq,PartialEq,Ord,PartialOrd,Hash,Serialize,Deserialize,Clone,Debug)]
pub struct BotCollisionHash(u32);

impl BotCollisionHash{
    fn new<T>(base:Fo<T>,a:&T,b:&T)->BotCollisionHash{                
        
        let a=((a as *const _ as usize - base.0 as usize )/core::mem::size_of::<T>()) as u32;
        let b=((b as *const _ as usize - base.0 as usize )/core::mem::size_of::<T>()) as u32;
        
        
        let [a,b]=if a>b{ //TODO important that else defaults to not moving?
            [b,a]
        }else{
            [a,b]
        };
        
    
        BotCollisionHash((a<<16) | b)
    }
}



#[derive(Eq,PartialEq,PartialOrd,Ord,Hash,Serialize,Deserialize,Debug,Copy,Clone)]
struct WallCollisionHash{
    a:u16,
    dir:grid::CardDir
}

impl WallCollisionHash{
    fn new<T>(base:Fo<T>,a:&T,dir:grid::CardDir)->WallCollisionHash{
        let a=((a as *const _ as usize - base.0 as usize) / core::mem::size_of::<T>()) as u16; 
        WallCollisionHash{a,dir}
    }
}


struct Fo<T>(*const T);

unsafe impl<T> Send for Fo<T>{}
unsafe impl<T> Sync for Fo<T>{}

impl<T> Copy for Fo<T> {}
impl<T> Clone for Fo<T> {
    fn clone(&self) -> Self {
        *self
    }
}


#[derive(Clone,Hash,Serialize,Deserialize,Debug)]
pub struct CollisionVelocitySolver<K>{
    last_bot_col:BTreeMap<BotCollisionHash,K>,
    last_wall_col:BTreeMap<WallCollisionHash,K>
}

use core::convert::TryFrom;
impl<B> CollisionVelocitySolver<B>{

    #[inline(always)]
    pub fn inner_try_into<A: TryFrom<B>>(self) -> Result<CollisionVelocitySolver<A>, A::Error> {
        let CollisionVelocitySolver{last_bot_col,last_wall_col}=self;
        
        
        //let mut last_bot_col_new:Vec<(BotCollisionHash,A)>=Vec::new();
        let mut last_bot_col_new=BTreeMap::new();
        for (a,b) in last_bot_col.into_iter(){
            last_bot_col_new.insert(a,TryFrom::try_from(b)?);
        }

        let mut last_wall_col_new=BTreeMap::new();
        for (a,b) in last_wall_col.into_iter(){
            last_wall_col_new.insert(a,TryFrom::try_from(b)?);
        }

        

        Ok(CollisionVelocitySolver{last_bot_col:last_bot_col_new,last_wall_col:last_wall_col_new})
        


    }
}

impl CollisionVelocitySolver<f32>{
    pub fn new()->CollisionVelocitySolver<f32>{
        CollisionVelocitySolver{last_bot_col:BTreeMap::new(),last_wall_col:BTreeMap::new()}
    }



    pub fn solve<A:Axis,T:Send+Sync+core::fmt::Debug>(
        &mut self,
        radius:f32,
        grid_viewport:&grid::GridViewPort,
        walls:&grid::Grid2D,
        tree:&mut dinotree_alg::collectable::CollectableDinoTree<A,NotNan<f32>,T>,
        pos_func:impl Fn(&T)->&Vec2<f32> + Send + Sync +Copy,
        vel_func:impl Fn(&mut T)->&mut Vec2<f32> + Send + Sync +Copy){
        
        let diameter=radius*2.0;
        let diameter2=diameter*diameter;
        let bias_factor=0.3;
        let allowed_penetration=radius*0.2;
        let num_iterations=5;
        let base=Fo(tree.get_bots().as_ptr());

        let mut collision_list={
            let ka3 = &self.last_bot_col;
            tree.collect_intersections_list_par(move |a:&mut T,b:&mut T|{
                let p2=*pos_func(b);
                let p1=*pos_func(a);
                let offset=p2-p1;
                let distance2=offset.magnitude2();
                
                if distance2<diameter2{
                    let distance=distance2.sqrt();
                    
                    let offset_normal=if distance<0.0001{
                        //Make a random direction
                        let rot=p1.y*2.0-p1.x;
                        vec2(rot.cos(),rot.sin())
                    }else{
                        offset/distance
                    };

                    let separation=(diameter-distance)/2.0;
                    let bias=-bias_factor*(1.0/num_iterations as f32)*( (-separation+allowed_penetration).min(0.0));
                    
                    let hash=BotCollisionHash::new(base,a,b);
                    //let impulse=20.0;
                    //let impulse=if true{
                    let impulse=if let Some(&impulse)=ka3.get(&hash){ //TODO inefficient to check if its none every time
                    
                        let k=offset_normal*impulse;
                        
                        *vel_func(a)-=k;
                        *vel_func(b)+=k;
                        impulse
                    }else{
                        0.0
                    };

                    Some((offset_normal,bias,impulse))
                }else{
                    None
                }
            })
        };

        //Package in one struct
        //so that there is no chance of mutating it twice
        #[derive(Debug)]
        struct WallCollision{
            collisions:[Option<(f32,Vec2<f32>,grid::CardDir,f32)>;2],
        }

        
        let mut wall_collisions={
            let ka3 = &self.last_wall_col;
            
            tree.collect_all_par(move |rect,a|{
                
                let arr=duckduckgeo::grid::collide::is_colliding(&walls,&grid_viewport,rect.as_ref(),radius);
                let create_collision=|bot:&mut T,dir:grid::CardDir,seperation:f32,offset_normal:Vec2<f32>|{
                    let bias=-bias_factor*(1.0/num_iterations as f32)*( (-seperation+allowed_penetration).min(0.0));

                    let impulse=if let Some(&impulse)=ka3.get(&WallCollisionHash::new(base,bot,dir)){ //TODO inefficient to check if its none every time
                        let k=offset_normal*impulse;
                        *vel_func(bot)+=k;
                        impulse
                    }else{
                        0.0
                    };
                    (bias,offset_normal,dir,impulse)
                };
                match arr[0]{
                    Some((seperation,dir,offset_normal))=>{
                        
                        let wall=match arr[1]{
                            Some((seperation,dir,offset_normal))=>{
                                let seperation=seperation*2.0f32.sqrt(); //Since we are pushing diagonally dont want to over push.
                                let first=Some(create_collision(a,dir,seperation,offset_normal));
                                let second=Some(create_collision(a,dir,seperation,offset_normal));
                                WallCollision{collisions:[first,second]}
                            },
                            None=>{
                                let first=Some(create_collision(a,dir,seperation,offset_normal));
                                WallCollision{collisions:[first,None]}
                            }
                        };
                        Some(wall)
                    },
                    None=>{
                        None
                    }
                }
            })
        };
        

        
        for _ in 0..num_iterations{

            
            
            collision_list.for_every_pair_mut_par(tree,|a,b,&mut (offset_normal,bias,ref mut acc)|{
                
                let vel=*vel_func(b)-*vel_func(a);

                let mass=0.2;
                let impulse=(bias-vel.dot(offset_normal))*mass;
                
                let p0=*acc;
                *acc=(p0+impulse).max(0.0);
                let impulse=*acc-p0;
                
                let k=offset_normal*impulse;
                //assert!(!k.is_nan());
                //println!("{:?}",k);
                *vel_func(a)-=k;
                *vel_func(b)+=k;
            });
                 
            
            wall_collisions.for_every_mut_par(tree,|bot,wall|{
                
                
                for k in wall.collisions.iter_mut(){
                    if let &mut Some((bias,offset_normal,_dir,ref mut acc))=k{
                        
                        let impulse=bias-vel_func(bot).dot(offset_normal);

                        let p0=*acc;
                        *acc=(p0+impulse).max(0.0);
                        let impulse=*acc-p0;

                        *vel_func(bot)+=offset_normal*impulse;
                    }
                }; 
            })
            
            
            
            
        } 
         


        self.last_bot_col.clear();
        self.last_wall_col.clear();

        let (ka2,ka3):(BTreeMap<_,_>,_)=rayon::join(||{
            
            collision_list.get(&tree).iter().flat_map(|a|a.iter()).map(|(a,b,(_,_,impulse))|{
                (BotCollisionHash::new(base,*a,*b),*impulse)
            }).collect()
            
        },
        ||{
            
            wall_collisions.get(&tree).iter().flat_map(|(bot,wall)|{
                let k=wall.collisions.iter().filter(|a|a.is_some()).map(|a|a.unwrap());
                k.map(move |(_,_,dir,impulse)|{
                    (WallCollisionHash::new(base,*bot,dir),impulse)
                })
            }).collect()
            
        });

        self.last_bot_col=ka2;
        self.last_wall_col=ka3;
    }
}
