# Game Component: Sprite

A basic sprite component for use in 2d games, with animations and directions.

## Installation

```bash
npm install @basementuniverse/sprite
```

## How to use

Create a sprite:

```ts
import { Sprite } from '@basementuniverse/sprite';

const sprite = new Sprite({
  // options here...
});
```

Update the sprite every frame:

```ts
sprite.update(dt);
```

Set the sprite's properties at any time:

```ts
sprite.position = vec(100, 100);
sprite.scale = 1;
sprite.rotation = 0; // radians
```

_(see [here](https://www.npmjs.com/package/@basementuniverse/vec) for `vec` library)_

Get or set the sprite's direction and animation:

```ts
sprite.animation = 'walk';
sprite.direction = 'right';
```

Start/pause/reset the current animation:

```ts
sprite.playAnimation();
sprite.pauseAnimation();
sprite.resetAnimation();
```

Render the sprite every frame:

```ts
sprite.draw(context);
```

## Options

```ts
export type SpriteOptions = {
  /**
   * The position of the sprite
   *
   * Defaults to (0, 0)
   */
  position?: vec;

  /**
   * The base size of the sprite
   *
   * If omitted, use the base image size, or fall back to the size of the
   * first image found in available directions for the default animation
   *
   * If a size still can't be found, default to (0, 0)
   */
  size?: vec;

  /**
   * Origin offset from top-left corner, used for rotation and scaling
   *
   * If omitted, this will be placed in the center of the sprite (based on
   * size, noting that size will be calculated first)
   *
   * @see SpriteOptions.size
   */
  origin?: vec;

  /**
   * The scale factor of the sprite
   *
   * Default is 1
   */
  scale?: number;

  /**
   * The sprite rotation, measured in radians
   *
   * Default is 0
   */
  rotation?: number;

  /**
   * An array of valid direction names
   *
   * By default a sprite will have one available direction: 'default'
   */
  directions: string[];

  /**
   * The initial direction of the sprite
   *
   * Default is 'default'
   */
  defaultDirection: string;

  /**
   * An optional base image
   *
   * If an animation frame doesn't exist (for example an animation has been
   * configured with 5 frames but only 3 images exist in the animation's
   * images array), we fall back to this image
   *
   * If an animation frame image can't be found and this base image doesn't
   * exist, we can still render the sprite - it will just have an empty image
   *
   * (this can be useful for sprites that act purely as hosts for attachment
   * points, for example)
   */
  image?: HTMLImageElement | HTMLCanvasElement;

  /**
   * A dictionary of animations
   *
   * Each animation can have multiple variants for different directions
   *
   * Use '*' as a direction name to indicate that a variant is available for
   * all directions and can be used as a fallback if the current direction
   * can't be found
   */
  animations: {
    [name: string]: {
      [direction: string]: SpriteAnimationOptions;
    };
  };

  /**
   * The initial animation for the sprite
   */
  defaultAnimation: string;

  /**
   * A list of attachment points
   */
  attachmentPoints?: SpriteAttachmentPointOptions[];

  /**
   * Optional hook called before rendering the sprite image
   */
  preRender?: (
    context: CanvasRenderingContext2D,
    sprite: Sprite
  ) => void;

  /**
   * Optional hook called after rendering the sprite image
   */
  postRender?: (
    context: CanvasRenderingContext2D,
    sprite: Sprite
  ) => void;

  /**
   * Optional debug options
   *
   * Can be a boolean value (in which case all sub-options will be set to the
   * same value), or an object allowing specific debug options to be enabled
   * individually
   */
  debug?: Partial<SpriteDebugOptions> | boolean;
};
```

_(see `build/index.d.ts` for more details)_

## Attachment points

Sprites can have multiple named attachment points. These attachment points can be given a default position, and the position can be modified via animation keyframes in each animation.

Fetch an attachment point by name:

```ts
sprite.getAttachmentPoint('weapon-left-hand');
// returns a vec: { x: number; y: number } or null
```

This could be useful for connecting multiple sprites together.
